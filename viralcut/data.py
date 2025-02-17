'''
data.py

This file acts as an interface data from the NCBI Datasets v2 REST API (dataset.py) and the local cache (cache.py).
The main file will make request and this interface will locate the data and perform some transformation.
'''

import json
import pickle
import os
import subprocess
import multiprocessing
import re
import zipfile as zf
from glob import glob
from io import BytesIO, TextIOWrapper

from . import config
from . import dataset
from . import cache
from . import utils
from .guide import Guide
from .collection import ViralCutCollection


def collection_from_pickle(filename):
    if config.VERBOSE:
        print(f'Loading ViralCutCollection from: {filename}')

    collection = ViralCutCollection()
    with open(filename, 'rb') as fp:
        collection = pickle.load(fp)

    if config.VERBOSE:
        print('Setting autosave path')
    collection._pickle_filepath = filename

    print('Loaded')
    return collection

def download_ncbi_genes(gene_ids):
    '''Downloads the sequence of the specified gene

    Arguments:
        gene_ids (list): List of NCBI gene IDs as integers

    Returns:
        A list of downloaded genes by ID
    '''

    # Convert gene_ids to list if non-list type is provided
    if not isinstance(gene_ids, list):
        gene_ids = [int(gene_ids)]

    # Some genes may already be cached,
    genes_to_download = cache.get_missing_genes(gene_ids)

    # which means we may not have anything to download.
    if not len(genes_to_download):
        if config.VERBOSE:
            print('No genes to download')
        return []

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(gene_ids) != len(genes_to_download):
        print((
            f'{len(gene_ids) - len(genes_to_download)} of {len(gene_ids)} requested exist in the '
            f'cache already.'
        ))

    # Keep track of which genes were actually downloaded.
    gene_ids_downloaded = []

    NCBI_gene_files = dataset.get_genes_by_id(genes_to_download)

    # Process the downloaded data without writing it to disk, yet.
    in_memory = BytesIO(NCBI_gene_files)
    with zf.ZipFile(in_memory) as zfp:
        files = zfp.namelist()

        if 'ncbi_dataset/data/gene.fna' not in files:
            raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/gene.fna`')

        # Extract the data report
        if 'ncbi_dataset/data/data_report.jsonl' not in files:
            raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/data_report.jsonl`')

        data_report_by_gene_id = {}

        with zfp.open('ncbi_dataset/data/data_report.jsonl') as fp:
            for line in fp:
                report = json.loads(line.strip())
                id = report['geneId']
                data_report_by_gene_id[id] = report

        # ZipFile.open() reads as a binary stream but we need a text stream.
        # https://docs.python.org/3/library/io.html#io.TextIOWrapper
        with TextIOWrapper(zfp.open('ncbi_dataset/data/gene.fna')) as fp:
            for header, seq in utils.parse_fna(fp):
                header_gene_id = get_properties_from_ncbi_fasta_header(header, key='GeneID')

                filename = f"{header_gene_id}.fna"
                cache_dir = cache.add_gene(header_gene_id)
                cached_file = os.path.join(cache_dir, filename)
                cache_data_report = os.path.join(cache_dir, 'data_report.json')

                with open(cached_file, 'w') as fp:
                    fp.write(f'{header}\n')
                    fp.write(f'{seq}\n')

                with open(cache_data_report, 'w') as fpW:
                    json.dump(data_report_by_gene_id[header_gene_id], fpW)

                if config.VERBOSE:
                    print(f'Downloaded to: {cache_dir}')

                gene_ids_downloaded.append(header_gene_id)

    return gene_ids_downloaded


def download_ncbi_assemblies(accessions, keep_exts=['fna'], merge=False):
    '''Download files associated with the provided accessions to a cache.

    Arguments:
        accessions (list): List of NCBI accessions as strings
        keep_exts (list):  Files with these extensions will be cached
        merge (bool):      Only cache files that have not been cached already

    Return:
        A list of accessions actually downloaded

    '''

    # Convert gene_ids to list if a string is provided
    if isinstance(accessions, str):
        accessions = [accessions]

    # Some accessions may already be cached,
    accs_to_download = cache.get_missing_assemblies(accessions)

    # which means we may not have anything to download.
    if not len(accs_to_download):
        if config.VERBOSE:
            print('No assemblies to download')
        return []

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(accessions) != len(accs_to_download):
        print((
            f'{len(accessions) - len(accs_to_download)} of {len(accessions)} '
            f'requested exist in the cache already.'
        ))

    # Keep track of which accessions were actually downloaded.
    accessions_downloaded = []

    for start in range(0, len(accs_to_download), config.NCBI_BATCH_SIZE):

        assembly_files = dataset.get_assembly_by_accession(accs_to_download[start:start+config.NCBI_BATCH_SIZE])
        
        # Process the downloaded data without writing it to disk, yet.
        in_memory = BytesIO(assembly_files)
        with zf.ZipFile(in_memory) as zfp:
            files = zfp.namelist()

            # Extract the data report
            if 'ncbi_dataset/data/assembly_data_report.jsonl' not in files:
                raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/assembly_data_report.jsonl`')

            data_report_by_accession = {}

            with zfp.open('ncbi_dataset/data/assembly_data_report.jsonl') as fp:
                for line in fp:
                    report = json.loads(line.strip())
                    accs = report['accession']
                    data_report_by_accession[accs] = report

            for path in files:

                # Be selective in which files get cached.
                accs = os.path.basename(os.path.dirname(path))
                filename = os.path.basename(path)
                cache_dir = cache.add_assembly(accs)
                cached_file = os.path.join(cache_dir, filename)
                cache_data_report = os.path.join(cache_dir, 'data_report.json')

                if not('GCF' in accs or 'GCA' in accs):
                    continue

                if filename.split('.')[-1] in keep_exts:

                    if not os.path.exists(cached_file):
                        with zfp.open(path) as fp, open(cached_file, 'wb') as fpW:
                            fpW.writelines(fp.readlines())

                        if config.VERBOSE:
                            print(f'Downloaded to: {cached_file}')

                    if not os.path.exists(cache_data_report):
                        with open(cache_data_report, 'w') as fpW:
                            json.dump(data_report_by_accession[accs], fpW)

                        if config.VERBOSE:
                            print(f'Downloaded to: {cache_data_report}')

                    accessions_downloaded.append(accs)

    # Done!
    return accessions_downloaded

def get_properties_from_ncbi_fasta_header(header, key=None):
    '''The NCBI Datasets API provides a multi-FASTA file when downloading genes.
    The headers are formatted like this:
        >NC_045512.2:26245-26472 E [organism=coronavirus 2] [gene_id=43740570] [chromosome=]
    This method will extract the value of the specified key from a FASTA header.

    In addition to the properties provided in square brackets, the following are also supported:
        - `accession`
        - `start`
        - `end`
        - `name`

    Arguments:
        header (string): The NCBI gene FASTA-header to parse
        key (string):    The key of the field to extract from the header.
                         If None, a dictionary of all properties is returned.

    Returns:
        See argument `key`.
    '''
    if header[0] == '>':
        header = header[1:]

    props = {}
    for i, prop in enumerate(header.split(' ')):
        if i == 0:
            # e.g., `NC_045512.2:26245-26472`
            acces, position = prop.split(':')
            start, end = position.split('-')
            props['accession'] = acces
            props['start'] = start
            props['end'] = end

        elif i == 1:
            # `E`
            props['name'] = prop

        else:
            # e.g., `[gene_id=43740570]`
            if prop[0] == '[' and prop[-1] == ']':
                k, v = prop[1:-1].split('=')
                props[k] = v

    if key is None:
        return props
    else:
        if key in props:
            return props[key]
        raise ValueError(f'Could not find `{key}` in header: `{header}`')


def _create_issl_index(commands):
    success = True
    for cmd in commands:
        try:
            subprocess.run(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            # subprocess.run(cmd)
        except Exception as e:
            if config.VERBOSE:
                print(f'Failed to run: {" ".join(cmd)}')
            success &= False
    return success

def create_issl_indexes(accessions,
    bin_extract=None,
    bin_issl=None,
    force=False,
    processors=0
):
    '''Extract offtargets and create ISSL index for each FNA file in the provided accession

    Arguments:
        accessions (list):     A list of accessions to process
        bin_extract (string):  A path to the extractOfftargets utility available in Crackling.
                               Set in configuration.
        bin_issl (string):     A path to the isslCreateIndex binary from Crackling.
                               Set in configuration.
        force (bool):          Rerun even if their output files already exist

    Returns:
        A list of accessions which an index was created for
    '''

    # Convert gene_ids to list if a string is provided
    if isinstance(accessions, str):
        accessions = [accessions]

    # Set defaults
    if bin_extract is None:
        bin_extract = config.BIN_EXTRACT

    if bin_issl is None:
        bin_issl = config.BIN_ISSL_IDX

    # remove accessions if the index already exists
    if not force:
        accessions = cache.get_missing_issl_index(accessions)

    args = []
    for accession in accessions:
        commands = []

        # For the .fna file associated with the accession
        fna_path = cache.get_assembly_fna(accession)
        fna = os.path.basename(fna_path)
        parent = os.path.pardir
        ots_file = os.path.join(parent, f"{fna}.offtargets.txt")
        issl_file = os.path.join(parent, f"{fna}.offtargets.issl")

        if not os.path.exists(ots_file) or force:
            # Extract off-targets
            commands.append([
                bin_extract,
                ots_file,
                fna_path
            ])

        if not os.path.exists(issl_file) or force:
            # Create ISSL index
            commands.append([
                bin_issl,
                ots_file,
                '20',
                '8',
                issl_file
            ])
        args.append(commands)

    errors = []
    if len(args) > 0:
        args, accessions = [(arg, acc) for arg, acc in zip(args, accessions) if len(arg) > 0]
        if processors == 1:
            for idx, arg in enumerate(args):
                if not _create_issl_index(arg):
                    errors.append(accessions[idx])
        else:
            with multiprocessing.Pool(os.cpu_count() if not processors else processors) as p:
                success = p.starmap(_create_issl_index, [[x] for x in args])
                errors = [accessions[idx] for idx, result in enumerate(success) if not result]
    return errors


def get_accession_from_tax_id(tax_ids):
    '''Given a list of taxonomy IDs, return a list of accessions.
    
    Arguments:
        tax_ids (list): A list of accessions
    
    Returns:
        A list of tax IDs, '-' indicates that an accession could not be found for tthe given taxon
    '''
    accession = ['-'] * len(tax_ids)

    for i in range(0, len(tax_ids), config.NCBI_BATCH_SIZE):
        print(f'{i}:{i+config.NCBI_BATCH_SIZE}')
    for i in range(0, len(tax_ids), config.NCBI_BATCH_SIZE):
        batch = tax_ids[i:i+config.NCBI_BATCH_SIZE]
        # Check GenBank first
        reports = dataset.get_genbank_dataset_reports_by_taxon(batch)
        for report in reports:
            accession[i+batch.index(report['organism']['tax_id'])] = report['accession']
            # accession[report['organism']['tax_id']] = report['accession']
        # Check Refseq for any entries not found in GenBank
        missed_tax_ids = [batch[idx] for idx, val in enumerate(batch) if val == '-']
        if len(missed_tax_ids) < 1:
            continue
        # missed_tax_ids = [t_id for t_id in batch if t_id not in accession.keys()]
        reports = dataset.get_refseq_dataset_reports_by_taxon(missed_tax_ids)
        for report in reports:
            accession[tax_ids.index(report['organism']['tax_id'])] = report['accession']
            # accession[report['organism']['tax_id']] = report['accession']

    missed_tax_ids = [tax_ids[idx] for idx, val in enumerate(accession) if val == '-']
    # missed_tax_ids = [t_id for t_id in tax_ids if t_id not in accession.keys()]
    if len(missed_tax_ids) > 0:
        print('Missed these boss')
        print(missed_tax_ids)

    return accession

def get_tax_ids_from_accessions(accessions, uniq=True):
    '''Given a list of accessions, return a list of taxonomy IDs.
    
    Arguments:
        accessions (list): A list of accessions
        uniq (bool):        If true, a set of tax IDs will be returned. Else, 
            duplicate tax IDs may be returned.
    
    Returns:
        A list or set of tax IDs (see arg `uniq`).
    '''
    tax_ids = []
    for accs in accessions:
        report = cache.get_assembly_report(accs)
        if 'taxId' in report['organism']:
            tax_ids.append(report['organism']['taxId'])
        elif config.VERBOSE:
            print(f'Could not find taxId of {accs}')
    return tax_ids

def get_guides_from_gene(gene_id):
    '''
    This method will extract all the guides from the gene
    
    '''
    pattern_forward = r'(?=([ATCG]{21}GG))'
    pattern_reverse = r'(?=(CC[ACGT]{21}))'

    guides = []
    seqs = cache.get_gene_seq(gene_id)
    for seq in seqs:
        for pattern, strand, seqModifier in [
            [pattern_forward, '+', lambda x : x],
            [pattern_reverse, '-', lambda x : utils.rc(x)]
        ]:
            p = re.compile(pattern)
            for m in p.finditer(seq):
                front_offset = 4 if strand == '+' else 3
                back_offset = 3 if strand == '+' else 4
                target30 = seqModifier(seq[m.start() - front_offset: m.start() + 23 + back_offset])
                target23 = target30[4:-3]
                guides.append((Guide(target23), target30, m.start(), strand))

    return guides


def get_guides_from_genome(accession):
    '''
    This method will extract all the guides from the gene
    
    '''
    pattern_forward = r'(?=([ATCG]{21}GG))'
    pattern_reverse = r'(?=(CC[ACGT]{21}))'

    guides = []
    seqs = cache.get_assembly_seq(accession)
    for seq in seqs:
        for pattern, strand, seqModifier in [
            [pattern_forward, '+', lambda x : x],
            [pattern_reverse, '-', lambda x : utils.rc(x)]
        ]:
            p = re.compile(pattern)
            for m in p.finditer(seq):
                front_offset = 4 if strand == '+' else 3
                back_offset = 3 if strand == '+' else 4
                target30 = seqModifier(seq[m.start() - front_offset: m.start() + 23 + back_offset])
                target23 = target30[4:-3]
                guides.append((Guide(target23), target30, m.start(), strand))

    return guides

def create_collection_from_gene(gene_id, guides):
    ''' this method will take a gene id and guides and create a new viral cut collection'''
    collection = ViralCutCollection()
    collection.target = ('gene', gene_id)
    collection.target_properties = cache.get_gene_report(gene_id)
    for guide, target30, start, strand in guides:
        if guide.seq not in collection:
            collection[guide.seq] = guide
            collection[guide.seq]['start'] = [start]
            collection[guide.seq]['end'] = [start + 23]
            collection[guide.seq]['strand'] = [strand]
            collection[guide.seq]['30mer'] = [target30]
            collection[guide.seq]['occurrences'] = 1
        else:
            collection[guide.seq]['start'].append(start)
            collection[guide.seq]['end'].append(start + 23)
            collection[guide.seq]['strand'].append(strand)
            collection[guide.seq]['30mer'].append(target30)
            collection[guide.seq]['occurrences'] += 1
    return collection

def create_collection_from_accession(accession, guides):
    ''' this method will take a accession and guides and create a new viral cut collection'''
    collection = ViralCutCollection()
    collection.target = ('accession', accession)
    collection.target_properties = cache.get_assembly_report(accession)
    for guide, target30, start, strand in guides:
        if guide.seq not in collection:
            collection[guide.seq] = guide
            collection[guide.seq]['start'] = [start]
            collection[guide.seq]['end'] = [start + 23]
            collection[guide.seq]['strand'] = [strand]
            collection[guide.seq]['30mer'] = [target30]
            collection[guide.seq]['occurrences'] = 1
        else:
            collection[guide.seq]['start'].append(start)
            collection[guide.seq]['end'].append(start + 23)
            collection[guide.seq]['strand'].append(strand)
            collection[guide.seq]['30mer'].append(target30)
            collection[guide.seq]['occurrences'] += 1
    return collection