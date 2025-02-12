import json
import pickle
import os
import subprocess
import multiprocessing
import zipfile as zf
from glob import glob
from io import BytesIO, TextIOWrapper

from . import config
from . import dataset
from . import cache
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

    ### # Recreate parts of the Collection that would not have been Pickled
    ### if config.VERBOSE:
    ###     print('Recreating components that were not Pickled')
    ### collection._ncbi = NCBITaxa()
    ###
    ### try:
    ###     if config.VERBOSE:
    ###         print('Trying to load tree from Newick')
    ###     collection._ncbi_tree = read_newick(collection._ncbi_tree_newick)
    ###
    ###     if config.VERBOSE:
    ###         print('Annotating tree')
    ###     collection._ncbi.annotate_tree(collection._ncbi_tree)
    ###
    ###     if config.VERBOSE:
    ###         print('Success')
    ###
    ### except Exception as e:
    ###     print(e)
    ###     print('... trying another way to load the tree')
    ###     collection._ncbi_tree = collection._ncbi.get_topology(
    ###         collection.get_tax_ids_in_analysis(),
    ###         intermediate_nodes=True
    ###     )
    ### #else:
    ### #    print('... trying even another way to load the tree')
    ### #    collection._ncbi_tree = collection._ncbi.get_topology(
    ### #        [config.ROOT_TAX_ID],
    ### #        intermediate_nodes=True
    ### #    )

    print('Loaded')
    return collection

def parse_fna(stream):
    '''Parse some iterable object as a multi-FASTA file.
    Yield after reading each FASTA block.

    Arguments:
        stream (iterable):  An iterable object to read
    '''
    header = None
    seqs = []
    for line in stream:
        line = line.strip()

        if line[0] == '>':
            if header is not None:
                yield header, ''.join(seqs)
            header = line
        else:
            seqs.append(line)
    yield header, ''.join(seqs)


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
            for header, seq in parse_fna(fp):
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



# def get_cached_genomes_by_id(gene_ids, download_missing=False):
#     '''A generator method that provides the sequence of the requested genes.

#     Arguments:
#         gene_ids (list):         A list of gene IDs to use in the generator
#         download_missing (bool): If a gene hasn't been cached, download it.

#     '''

#     if download_missing:
#         to_download = []

#         for gene_id in gene_ids:
#             if not os.path.exists(get_gene_cache(gene_id)):
#                 to_download.append(gene_id)

#         download_ncbi_genes(to_download)

#     for gene_id in gene_ids:
#         fna_fp = os.path.join(get_gene_cache(gene_id), f"{gene_id}.fna")

#         with open(fna_fp, 'r') as fp:
#             for header, seq in parse_fna(fp):
#                 yield gene_id, seq


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
        accessions = cache.get_missing_issl_index()

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

    args, accessions = [arg, acc for arg, acc in zip(args, accessions) if len(arg) > 0]
    errors = []
    if processors == 1:
        for idx, arg in enumerate(args):
            if not _create_issl_index(arg):
                errors.append(accessions[idx])
    else:
        with multiprocessing.Pool(os.cpu_count() if not processors else processors) as p:
            success = p.starmap(_create_issl_index, [[x] for x in args])
            errors = [accessions[idx] for idx, result in enumerate(success) if not result]:
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
