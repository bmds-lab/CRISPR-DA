'''
data.py

This file acts as an interface data from the NCBI Datasets v2 REST API (dataset.py) and the local cache (cache.py).
The main file will make request and this interface will locate the data and perform some transformation.
'''
import re
import os
import glob
import json
import heapq
import shutil
import pickle
import tempfile
import subprocess
import multiprocessing
import zipfile as zf
from pathlib import Path
from importlib import resources
from io import BytesIO, TextIOWrapper
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa

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
    genes_to_download = cache.get_missing_files(gene_ids, '.fna')

    # which means we may not have anything to download.
    if not len(genes_to_download):
        if config.VERBOSE:
            print('No genes to download')
        return gene_ids

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(gene_ids) != len(genes_to_download):
        print((
            f'{len(gene_ids) - len(genes_to_download)} of {len(gene_ids)} requested exist in the '
            f'cache already.'
        ))

    # Keep track of which genes were actually downloaded.
    gene_ids_downloaded = [x for x in gene_ids if x not in genes_to_download]

    success, NCBI_gene_files = dataset.get_genes_by_id(genes_to_download)
    if not success:
        raise RuntimeError('Unable to download requested genes')

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

                cache_dir = cache.add_entry(header_gene_id)
                cached_file = cache_dir / f'{header_gene_id}.fna'
                cache_data_report = cache_dir / 'data_report.json'

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
    accs_to_download = cache.get_missing_files(accessions, '.fna')

    # which means we may not have anything to download.
    if not len(accs_to_download):
        if config.VERBOSE:
            print('No assemblies to download')
        return accessions

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(accessions) != len(accs_to_download):
        print((
            f'{len(accessions) - len(accs_to_download)} of {len(accessions)} '
            f'requested exist in the cache already.'
        ))

    # Keep track of which accessions were actually downloaded.
    accessions_downloaded = [x for x in accessions if x not in accs_to_download]
    for start in range(0, len(accs_to_download), config.NCBI_BATCH_SIZE):

        success, assembly_files = dataset.get_assembly_by_accession(accs_to_download[start:start+config.NCBI_BATCH_SIZE])
        if not success:
            raise RuntimeError('Unable to download requested assemblies')
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

            for file in [Path(x) for x in files]:

                # Be selective in which files get cached.
                accs = file.parent.name
                filename = file.name
                if not('GCF' in accs or 'GCA' in accs):
                    continue
                cache_dir = cache.add_entry(accs)
                cached_file = cache_dir / filename
                cache_data_report = cache_dir / 'data_report.json'

                if filename.split('.')[-1] in keep_exts:

                    if not cached_file.exists():
                        with zfp.open(str(file)) as fp, open(cached_file, 'wb') as fpW:
                            fpW.writelines(fp.readlines())

                        if config.VERBOSE:
                            print(f'Downloaded to: {cached_file}')

                    if not cache_data_report.exists():
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

def extract_offtargets(input, output, max_open_files=1000):
    try:
        pattern_forward_offsite = r"(?=([ACGT]{21}[AG]G))"
        pattern_reverse_offsite = r"(?=(C[CT][ACGT]{21}))"
        # Create multiprocessing pool
        # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.starmap
        tempDir = tempfile.TemporaryDirectory()
        with open(str(input), 'r') as inFile:
            for header, seq in utils.parse_fna(inFile):
                with tempfile.NamedTemporaryFile('w', delete=False, dir=tempDir.name, suffix='_split') as outFile:
                    outFile.write(f'{header}\n')
                    outFile.write(f'{seq}\n')

        explodedFiles = [file for file in glob.glob(f'{tempDir.name}/*_split')]
        for explodedFile in explodedFiles:
            with open(explodedFile, 'r') as inFile:
                lines = inFile.readlines()
            if len(lines) != 2:
                raise RuntimeError("Error: expected the file to have two lines (>header,seq)")
            header = lines[0].strip()
            seq = lines[1].strip()
            with tempfile.NamedTemporaryFile('w', delete=False, dir=tempDir.name, suffix='_unsorted') as outFile:
                for strand, pattern, seqModifier in [
                    ['positive', pattern_forward_offsite, lambda x : x],
                    ['negative', pattern_reverse_offsite, lambda x : utils.rc(x)]
                ]:
                    match_chr = re.findall(pattern, seq)
                    for i in range(0,len(match_chr)):
                        outFile.write(f'{seqModifier(match_chr[i])[0:20]}\n')

        for unsortedFile in glob.glob(f'{tempDir.name}/*_unsorted'):
            with tempfile.NamedTemporaryFile('w', delete=False, dir=tempDir.name, suffix='_sorted') as outFile, open(unsortedFile, 'r') as inFile:
                lines = inFile.readlines()
                lines.sort()
                outFile.writelines(lines)

        sortedFiles = [file for file in glob.glob(f'{tempDir.name}/*_sorted')]
        while len(sortedFiles) > 1:
            mergedFile = tempfile.NamedTemporaryFile(delete = False)
            while True:
                try:
                    sortedFilesPointers = [open(file, 'r') for file in sortedFiles[:max_open_files]]
                    break
                except OSError as e:
                    if e.errno == 24:
                        utils.printer(f'Attempted to open too many files at once (OSError errno 24)')
                        max_open_files = max(1, int(max_open_files / 2))
                        utils.printer(f'Reducing the number of files that can be opened by half to {max_open_files}')
                        continue
                    raise e
            with open(mergedFile.name, 'w') as f:
                f.writelines(heapq.merge(*sortedFilesPointers))
            for file in sortedFilesPointers:
                file.close()
            sortedFiles = sortedFiles[max_open_files:] + [mergedFile.name]

        shutil.move(sortedFiles[0], output)
        tempDir.cleanup()
        return True
    except Exception as e:
        if config.VERBOSE:
            print(f'Failed to extract off-targets from: {input}')
            print(e)
            return False

def create_issl_indexes(accessions, force=True, processors=os.cpu_count()):
    '''Extract offtargets and create ISSL index for each FNA file in the provided accession

    Arguments:
        accessions (list):     A list of accessions to process
        force (bool):          Rerun even if their output files already exist

    Returns:
        A list of accessions which an index was created for
    '''

    # Convert gene_ids to list if a string is provided
    if isinstance(accessions, str):
        accessions = [accessions]

    with resources.path('viralcut.resources', 'slice-4-5.txt') as resource:
        sliceFile = resource

    with resources.path('viralcut.resources', 'ISSLCreateIndex') as resource:
        isslCreateIndexBin = resource

    # remove accessions if the index already exists
    if not force:
        accessions = cache.get_missing_files(accessions, '.issl')

    args = []
    for accession in accessions:
        fnaFile = cache.get_file(accession, '.fna')
        offtargetFile = fnaFile.parent / f"{fnaFile.stem}_offtargets.txt"
        if (not offtargetFile.exists()) or force:
            # Extract off-targets
            args.append([
                fnaFile,
                offtargetFile,
                100
            ])

    extractOfftargetErrors = []
    args = [(arg, acc) for arg, acc in zip(args, accessions) if len(arg) > 0]
    if len(args) > 0:
        args, accessions = zip(*args)
        if processors == 1:
            for idx, arg in enumerate(args):
                if not extract_offtargets(*arg):
                    extractOfftargetErrors.append(accessions[idx])
        else:
            with multiprocessing.Pool(processors) as p:
                success = p.starmap(extract_offtargets, args)
                extractOfftargetErrors = [accessions[idx] for idx, result in enumerate(success) if not result]

    args = []
    for accession in accessions:
        fnaFile = cache.get_file(accession, '.fna')
        offtargetFile = fnaFile.parent / f"{fnaFile.stem}_offtargets.txt"
        isslIndexFile = fnaFile.parent / f"{fnaFile.stem}.issl"
        if (not isslIndexFile.exists()) or force:
            # Create ISSL index
            args.append([
                isslCreateIndexBin,
                offtargetFile,
                sliceFile,
                '20',
                isslIndexFile
            ])

    buildISSLIndexErrors = []
    args = [(arg, acc) for arg, acc in zip(args, accessions) if len(arg) > 0]
    if len(args) > 0:
        args, accessions = zip(*args)
        if processors == 1:
            for idx, arg in enumerate(args):
                if not utils.run_command(arg):
                    buildISSLIndexErrors.append(accessions[idx])
        else:
            with multiprocessing.Pool(processors) as p:
                success = p.starmap(utils.run_command, [[x] for x in args])
                buildISSLIndexErrors = [accessions[idx] for idx, result in enumerate(success) if not result]
    return extractOfftargetErrors + buildISSLIndexErrors


def get_accession_from_tax_id(tax_ids):
    '''Given a list of taxonomy IDs, return a list of accessions.
    
    Arguments:
        tax_ids (list): A list of accessions
    
    Returns:
        A list of tax IDs, '-' indicates that an accession could not be found for the given taxon
    '''
    accession = ['-'] * len(tax_ids)

    for i in range(0, len(tax_ids), config.NCBI_BATCH_SIZE):
        batch = tax_ids[i:i+config.NCBI_BATCH_SIZE]
        # Check GenBank first
        success, reports = dataset.get_genbank_dataset_reports_by_taxon(batch)
        if not success:
            raise RuntimeError('Unable to download requested assemblies reports')
        for report in reports:
            accession[i+batch.index(report['organism']['tax_id'])] = report['accession']
        # Check Refseq for any entries not found in GenBank
        missed_tax_ids = [batch[idx] for idx, val in enumerate(batch) if val == '-']
        if len(missed_tax_ids) < 1:
            continue
        success, reports = dataset.get_refseq_dataset_reports_by_taxon(missed_tax_ids)
        if not success:
            raise RuntimeError('Unable to download requested assemblies reports')
        for report in reports:
            accession[tax_ids.index(report['organism']['tax_id'])] = report['accession']

    missed_tax_ids = [tax_ids[idx] for idx, val in enumerate(accession) if val == '-']

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
        reportFile = cache.get_file(accs, 'data_report.json')
        with open(reportFile, 'r') as inFile:
            report = json.load(inFile)
        if 'taxId' in report['organism']:
            tax_ids.append(report['organism']['taxId'])
        elif config.VERBOSE:
            print(f'Could not find taxId of {accs}')
    return tax_ids

def get_name_from_accessions(accessions):
    '''Given a list of accessions, return a list of scientific names.
    
    Arguments:
        accessions (list): A list of accessions
    
    Returns:
        A list of names (str).
    '''
    names = []
    for accs in accessions:
        reportFile = cache.get_file(accs, 'data_report.json')
        with open(reportFile, 'r') as inFile:
            report = json.load(inFile)
        if 'taxId' in report['organism']:
            names.append(report['organism']['organismName'])
        elif config.VERBOSE:
            print(f'Could not find taxId of {accs}')
    return names

def extract_guides(id):
    '''
    This method will extract all the guides from the gene or accession.

    Arguments:
        id (str): This either the gene id or accession number.
    
    Returns:
        guides list[tuple]: A list of extracted guides containing
                            A guide object, header, extended 30mer seq, start location and strand
    '''
    pattern_forward = r'(?=([ATCG]{25}GG[ATCG]{3}))'
    pattern_reverse = r'(?=([ATCG]{3}CC[ACGT]{25}))'

    guides = []
    with open(cache.get_file(id, '.fna'), 'r') as inFile:
        for header, seq in utils.parse_fna(inFile):
            for pattern, strand, seqModifier in [
                [pattern_forward, '+', lambda x : x],
                [pattern_reverse, '-', lambda x : utils.rc(x)]
            ]:
                p = re.compile(pattern)
                for m in p.finditer(seq):
                    target30 = seqModifier(seq[m.start(): m.start() + 30])
                    target23 = target30[4:-3]
                    guides.append((Guide(target23), header, target30, m.start(), strand))
    return guides

def create_collection(id, guides):
    ''' this method will take a gene id and guides and create a new viral cut collection'''
    collection = ViralCutCollection()
    collection.target = cache.get_file(id, 'data_report.json')
    for guide, header, target30, start, strand in guides:
        if guide.seq not in collection:
            collection[guide.seq] = guide
            collection[guide.seq]['header'] = [header]
            collection[guide.seq]['start'] = [start]
            collection[guide.seq]['end'] = [start + 23]
            collection[guide.seq]['strand'] = [strand]
            collection[guide.seq]['30mer'] = [target30]
            collection[guide.seq]['occurrences'] = 1
        else:
            collection[guide.seq]['header'].append(header)
            collection[guide.seq]['start'].append(start)
            collection[guide.seq]['end'].append(start + 23)
            collection[guide.seq]['strand'].append(strand)
            collection[guide.seq]['30mer'].append(target30)
            collection[guide.seq]['occurrences'] += 1
    return collection

def get_accession_from_root(root_tax_id):
        ncbi = NCBITaxa()
        tree = ncbi.get_topology([root_tax_id], intermediate_nodes=True)
        tax_ids = [node.taxid for node in tree.traverse("postorder")]
        accession = get_accession_from_tax_id(tax_ids)
        return tax_ids, [x for x in accession if x != '-']