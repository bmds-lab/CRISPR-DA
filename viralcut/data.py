import csv
import json
import os
import subprocess
import multiprocessing
import zipfile as zf

from io import BytesIO, TextIOWrapper

import ncbi.datasets
from glob import glob

from . import config


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

    # Convert gene_ids to list if a string is provided
    if isinstance(gene_ids, str):
        gene_ids = [gene_ids]

    # Some genes may already be cached,
    genes_to_download = [
        gene_id
        for gene_id in gene_ids
        if not os.path.exists(get_gene_cache(gene_id, mkdir=False))
    ]

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

    # New instance of the NCBI API interface.
    api_instance = ncbi.datasets.GeneApi(ncbi.datasets.ApiClient())

    api_response = api_instance.download_gene_package(
        [int(x) for x in gene_ids],
        include_annotation_type=['FASTA_GENE'],
        _preload_content=False
    )

    # Process the downloaded data without writing it to disk, yet.
    in_memory = BytesIO(api_response.data)
    with zf.ZipFile(in_memory) as zfp:
        files = zfp.namelist()

        if 'ncbi_dataset/data/gene.fna' not in files:
            raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/gene.fna`')

        # ZipFile.open() reads as a binary stream but we need a text stream.
        # https://docs.python.org/3/library/io.html#io.TextIOWrapper
        with TextIOWrapper(zfp.open('ncbi_dataset/data/gene.fna')) as fp:
            for header, seq in parse_fna(fp):
                headergene_id = get_properties_from_ncbi_fasta_header(header, key='GeneID')

                filename = f"{headergene_id}.fna"
                cache_dir = get_gene_cache(headergene_id)
                cached_file = os.path.join(cache_dir, filename)

                with open(cached_file, 'w') as fp:
                    fp.write(f'{header}\n')
                    fp.write(f'{seq}\n')

                gene_ids_downloaded.append(headergene_id)

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
    accs_to_download = [
        accs
        for accs in accessions
        if not os.path.exists(get_assembly_cache(accs, mkdir=False)) or merge
    ]

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

    # New instance of the NCBI API interface.
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

    for start in range(0, len(accs_to_download), config.NCBI_ACCESSION_BATCH_SIZE):
        api_response = api_instance.download_assembly_package(
            accs_to_download[start:start+config.NCBI_ACCESSION_BATCH_SIZE],
            _preload_content=False
        )

        # Process the downloaded data without writing it to disk, yet.
        in_memory = BytesIO(api_response.data)
        with zf.ZipFile(in_memory) as zfp:
            files = zfp.namelist()
            
            # Extract the data report
            if 'ncbi_dataset/data/assembly_data_report.jsonl' not in files:
                raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/assembly_data_report.jsonl`')

            dataReportByAccession = {}
            
            with zfp.open('ncbi_dataset/data/assembly_data_report.jsonl') as fp:
                for line in fp:
                    report = json.loads(line.strip())
                    accs = report['assemblyInfo']['assemblyAccession']
                    dataReportByAccession[accs] = report

            for path in files:

                # Be selective in which files get cached.
                accs = os.path.basename(os.path.dirname(path))
                filename = os.path.basename(path)
                cache_dir = get_assembly_cache(accs)
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
                            json.dump(dataReportByAccession[accs], fpW)
                        
                        if config.VERBOSE:
                            print(f'Downloaded to: {cache_data_report}')

                    accessions_downloaded.append(accs)

    # Done!
    return accessions_downloaded


def get_cached_gene_seqs_by_id(gene_ids, download_missing=False):
    '''A generator method that provides the sequence of the requested genes.

    Arguments:
        gene_ids (list):         A list of gene IDs to use in the generator
        download_missing (bool): If a gene hasn't been cached, download it.

    '''

    if download_missing:
        to_download = []

        for gene_id in gene_ids:
            if not os.path.exists(get_gene_cache(gene_id)):
                to_download.append(gene_id)

        download_ncbi_genes(to_download)

    for gene_id in gene_ids:
        fna_fp = os.path.join(get_gene_cache(gene_id), f"{gene_id}.fna")

        with open(fna_fp, 'r') as fp:
            for header, seq in parse_fna(fp):
                yield gene_id, seq


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
            # `NC_045512.2:26245-26472`
            acces, position = prop.split(':')
            start, end = position.split('-')
            props['accession'] = acces
            props['start'] = int(start)
            props['end'] = int(end)

        elif i == 1:
            # `E`
            props['name'] = prop

        else:
            # `[gene_id=43740570]`
            if prop[0] == '[' and prop[-1] == ']':
                k, v = prop[1:-1].split('=')
                props[k] = v

    if key is None:
        return props
    else:
        if key in props:
            return props[key]
        raise ValueError(f'Could not find `{key}` in header: `{header}`')


def get_assembly_cache(accession, mkdir=True):
    '''Generates a directory to store data for a particular accession.
    This function should be changed to your needs but keep in mind that changes
    may break the logic structure of your filesystem.

    Arguments:
        accession (string): The NCBI accession used for creating the directory
        mkdir (bool):       True, the directory is created.
                            False, the directory is not created.
                            Set to False if you just want to retrieve the path to the cache.

    Returns:
        The path to the directory where the NCBI data will be cached
    '''
    path = os.path.join(config.CACHE, accession[:config.CACHE_PREFIX_LENGTH], accession)
    if mkdir:
        os.makedirs(path, exist_ok=True)
    return path


def get_gene_cache(gene_id, mkdir=True):
    '''Generates a directory to store data for a particular gene.
    This function should be changed to your needs but keep in mind that changes
    may break the logic structure of your filesystem.

    Arguments:
        gene_id (string): The NCBI gene ID used for creating the directory
        mkdir (bool):     True, the directory is created.
                          False, the directory is not created.
                          Set to False if you just want to retrieve the path to the cache.

    Returns:
        The path to the directory where the NCBI data will be cached
    '''
    path = os.path.join(config.CACHE, f"gene-{gene_id}")
    if mkdir:
        os.makedirs(path, exist_ok=True)
    return path


def create_issl_indexes(accessions,
    bin_extract=None,
    bin_issl=None,
    stdout_extract=subprocess.DEVNULL,
    stdout_issl=subprocess.DEVNULL,
    stderr_extract=subprocess.DEVNULL,
    stderr_issl=subprocess.DEVNULL,
    force=False
):
    '''Extract offtargets and create ISSL index for each FNA file in the provided accession

    Arguments:
        accessions (list):     A list of accessions to process
        bin_extract (string):  A path to the extractOfftargets utility available in Crackling.
                               Set in configuration.
        bin_issl (string):     A path to the isslCreateIndex binary from Crackling.
                               Set in configuration.
        stdout_extract (file): File-like ojbect to redirect stdout for extractOfftargets to
        stdout_issl (file):    File-like ojbect to redirect stdout for isslCreateIndex to
        stderr_extract (file): File-like ojbect to redirect stderr for extractOfftargets to
        stderr_issl (file):    File-like ojbect to redirect stderr for isslCreateIndex to
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
    accessions = [
        accs
        for accs in accessions
        if len(
            glob(
                os.path.join(
                    get_assembly_cache(accs), 
                    '*.issl'
                )
            )
        ) == 0
    ]

    indexesCreateFor = []

    for accession in accessions:

        # Determine cache location
        cache = get_assembly_cache(accession)

        # Have all indexes been built successfully?
        success = True

        # Logic check: there should only be one *.fna per accession
        num_fna_files = len(glob(os.path.join(cache, '*.fna')))
        if num_fna_files != 1:
            raise RuntimeError((
                f'Found {num_fna_files} FNA files for {accession}. '
                f'There should be one.'
            ))

        # For the .fna file associated with the accession
        fna_path = glob(os.path.join(cache, '*.fna'))[0]

        fna = os.path.basename(fna_path)

        ots_file = os.path.join(cache, f"{fna}.offtargets.txt")
        issl_file = os.path.join(cache, f"{fna}.offtargets.issl")

        if not os.path.exists(ots_file) or force:

            # Extract off-targets
            cmd = [
                bin_extract,
                ots_file,
                fna_path
            ]

            try:
                subprocess.run(cmd,
                    stdout=stdout_extract,
                    stderr=stderr_extract
                )
            except Exception as e:
                if config.VERBOSE:
                    print(f'Failed to run: {" ".join(cmd)}')
                success &= False

        if not os.path.exists(issl_file) or force:
            # Create ISSL index
            cmd = [
                bin_issl,
                ots_file,
                '20',
                '8',
                issl_file
            ]
            try:
                subprocess.run(cmd,
                    stdout=stdout_issl,
                    stderr=stderr_issl
                )
            except Exception as e:
                if config.VERBOSE:
                    print(f'Failed to run: {" ".join(cmd)}')
                success &= False
        
        if success:
            indexesCreateFor.append(accession)
            
    return indexesCreateFor


def get_accessions_from_ncbi_table_export(filePath, **kwargs):
    '''The 'NCBI Genome Information by Organism' table is a good way to obtain a filtered list of
    NCBI accessions. This method parses the CSV file that is provided when clicking 'Download' on
    this page https://www.ncbi.nlm.nih.gov/genome/browse/. Please ensure the 'Assembly' column is
    present before exporting.
    
    Arguments:
        filePath (string): The path to the CSV file provided when Download is clicked on the webpage
        **kwargs: Keyword arguments forwarded to `csv.reader()`
        
    Returns:
        A list of NCBI accessions found in the file
    
    '''
    
    accessions = []
    
    # Find the `Assembly` column then extract all values in it
    with open(filePath, 'r') as fp:
        fp_rdr = csv.reader(fp, **kwargs)
        idx_assembly = None
        
        for i, row in enumerate(fp_rdr):
            if i == 0:
                if 'Assembly' in row:
                    idx_assembly = row.index('Assembly')
                    continue
                else:
                    raise ValueError('The column `Assembly` must be present in the file, yet it could not be found.')
            
            accessions.append(row[idx_assembly])
    
    return accessions