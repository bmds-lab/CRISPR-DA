
import os
import csv
import json
import shutil
from pathlib import Path
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


def get_gene_seq(gene_id):
    '''
    This function will get the gene seqeuence from the cache.

    Arguments:
        gene_id :   The gene ID to retrieve from the cache 

    Returns:
        Seqeunces : A list of sequences from the gene fna file
    '''
    gene_fna = Path(config.CACHE) / f"gene-{gene_id}" / f"{gene_id}.fna"
    seqs = []
    with open(gene_fna, 'r') as inFile: 
        for header, seq in parse_fna(inFile):
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

def add_assembly(accession):
    '''
    This function creates a new directory in the cache.
    It WILL remove the existing entry.
    
    Arguments:
        gene_id (int): NCBI gene IDs as integers
    '''
    cache = cache / accession[:config.CACHE_PREFIX_LENGTH] / accession
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir()

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

def get_cached_gene_information_by_id(gene_id):
    '''Given some gene identifier, collect as properties about the gene.
    Raises an exception if the gene does not exist in the cache.

    Arguments:
        gene_id (int):  The NCBI gene identifier

    Returns:
        A dictionary of properties
    '''
    dir = get_gene_cache(gene_id, mkdir=False)
    report_path = os.path.join(dir, 'data_report.json')
    if os.path.exists(report_path):
        with open(report_path, 'r') as fp:
            return json.loads(fp.readline())

    if config.VERBOSE:
        print(f'Could not find data report for gene #{gene_id}')

    return {}

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