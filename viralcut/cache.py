
import os
import csv
import json
import shutil
from pathlib import Path
from . import config
from . import utils


def get_missing_genes(gene_ids):
    '''
    This function will check if the provided gene ids exist in the cache.
    
    Arguments:
        gene_ids (list): List of NCBI gene IDs as integers

    Returns:
        A list of gene ids NOT downloaded in the cache.
    '''
    cache = Path(config.CACHE)
    missing_genes = []
    for gene_id in gene_ids:
        if not (cache / f"gene-{gene_id}").exists():
            missing_genes.append(gene_id)
    return missing_genes

def add_gene(gene_id):
    '''
    This function creates a new directory in the cache.
    It WILL remove the existing entry.
    
    Arguments:
        gene_id (int): NCBI gene IDs as integers
    '''
    cache = Path(config.CACHE) / f"gene-{gene_id}"
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir()
    return cache


def get_gene_seq(gene_id):
    '''
    This function will get the gene seqeuence from the cache.

    Arguments:
        gene_id :   The gene ID to retrieve from the cache 

    Returns:
        Seqeunces : A list of sequences from the gene fna file
    '''
    fna = Path(config.CACHE) / f"gene-{gene_id}" / f"{gene_id}.fna"
    seqs = []
    with open(fna, 'r') as inFile: 
        for header, seq in utils.parse_fna(inFile):
            seqs.append(seq)
    return seqs

def get_assembly_seq(accession: str):
    '''
    This function will get the gene seqeuence from the cache.

    Arguments:
        gene_id :   The gene ID to retrieve from the cache 

    Returns:
        Seqeunces : A list of sequences from the gene fna file
    '''
    fna = get_assembly_fna(accession)
    seqs = []
    with open(fna, 'r') as inFile: 
        for header, seq in utils.parse_fna(inFile):
            seqs.append(seq)
    return seqs

def get_missing_assemblies(accessions):
    '''
    This function will check if the provided accessions exist in the cache.
    
    Arguments:
        accessions (list): List of NCBI assembly accessions

    Returns:
        A list of accessions NOT downloaded in the cache.
    '''
    cache = Path(config.CACHE)
    missing_accessions = []
    for accession in accessions:
        if not (cache / accession[:config.CACHE_PREFIX_LENGTH] / accession).exists():
            missing_accessions.append(accession)
    return missing_accessions

def add_assembly(accession: str):
    '''
    This function creates a new directory in the cache.
    It WILL remove the existing entry.
    
    Arguments:
        gene_id (int): NCBI gene IDs as integers
    '''
    cache = Path(config.CACHE) / accession[:config.CACHE_PREFIX_LENGTH] / accession
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir()

def get_assembly_report(accession: str):
    '''
    Retrieves the data_report.json file for the given accesssion'''
    report = Path(config.CACHE) / accession[:config.CACHE_PREFIX_LENGTH] / accession / 'data_report.json'
    if not report.exists():
        raise RuntimeError(f'Could not find file {report} in cache')
    with open(report, 'r') as fp:
        data = json.load(fp)
    return data

def get_gene_report(gene_id: str):
    '''
    Retrieves the data_report.json file for the given accesssion'''
    report = Path(config.CACHE) / f'gene-{gene_id}' /'data_report.json'
    if not report.exists():
        raise RuntimeError(f'Could not find file {report} in cache')
    with open(report, 'r') as fp:
        data = json.load(fp)
    return data

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


def get_accession_issl_index(accession: str):
    '''
    This function will get the issl for the given accession
    
    Arguments:
        accessions (list): List of NCBI assembly accessions

    Returns:
        A list of accessions that do NOT have an ISSL index built in the cache
    '''

    cache = Path(config.CACHE) / accession[:config.CACHE_PREFIX_LENGTH] / accession
    files = [x for x in cache.glob('*.issl')]
    if len(files) == 0:
        raise RuntimeError(f'Could not find index for {accession}')
    return files[0]

def get_missing_issl_index(accessions: list[str]):
    '''
    This function will check if the provided accessions have an ISSL index in the cache
    
    Arguments:
        accessions (list): List of NCBI assembly accessions

    Returns:
        A list of accessions that do NOT have an ISSL index built in the cache
    '''
    missing_indexes = []
    for accession in accessions:
        cache = Path(config.CACHE) / accession[:config.CACHE_PREFIX_LENGTH] / accession
        if len([cache.glob('*.issl')]) == 0:
            missing_indexes.append(accession)
    return missing_indexes

def get_assembly_fna(accession: str):
    '''
    This function will get the associated FNA file for the given accession
    
    Arguments:
        accession : A NCBI assembly accession

    Returns:
        The path of the FNA file for the provide accession
    '''
    cache = Path(config.CACHE) / accession[:config.CACHE_PREFIX_LENGTH] / accession
    fna_files = [x for x in cache.glob('*.fna')]
    if len(fna_files) > 1:
            raise RuntimeError((
                f'Found {len(cache.glob('*.issl'))} FNA files for {accession}. '
                f'There should be one.'
            ))
    elif len(fna_files) == 0:
            raise RuntimeError((
                f'Could not find FNA files for {accession}. '
                f'There should be one.'
            ))
    return fna_files[0]

def get_accessions_from_ncbi_table_export(
    filePath=config.NCBI_HUMAN_VIRUSES_TABLE,
    **kwargs
):
    '''The 'NCBI Genome Information by Organism' table is a good way to obtain a filtered list of
    NCBI accessions. This method parses the CSV file that is provided when clicking 'Download' on
    this page https://www.ncbi.nlm.nih.gov/genome/browse/. Please ensure the 'Assembly' column is
    present before exporting.

    Arguments:
        filePath (string):  (optional) The path to the CSV file provided when Download is clicked on
                            the webpage. By default, the precompiled list provided with ViralCut
                            will be used.
        **kwargs:           Keyword arguments forwarded to `csv.reader()`

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

def get_cached_tax_ids():
    '''This function crawls the assemblies cache to extract taxonomy IDs

    Returns:
        A list of taxonomy IDs
    '''

    tax_ids = []

    for path, subdirs, files in os.walk(config.CACHE):
        for name in files:
            filename = os.path.basename(name)
            if filename == 'data_report.json':
                with open(os.path.join(path, filename), 'r') as fp:
                    report = json.loads(fp.readline())
                if 'taxId' in report:
                    tax_ids.append(report['taxId'])

    return tax_ids