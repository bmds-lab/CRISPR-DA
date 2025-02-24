'''
dataset.py

This file gets data from NCBI datasets using NCBI Datasets v2 REST API.
More details can be found here:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/

'''

import os
import netrc
from time import sleep
from requests import get
from requests.auth import HTTPBasicAuth

GENOME_ENDPOINT = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome'
GENE_ENDPOINT = 'https://api.ncbi.nlm.nih.gov/datasets/v2/gene'
AUTH = None

# Check for an NCBI API key saved as an enviroment variable or saved in the users .netrc file
# NOTE: This is not required but will help avoid request rate limits
try:
    # Attempt to get API key from enviroment variable 'NCBI_API_KEY'
    AUTH = HTTPBasicAuth('api-key', os.environ['NCBI_API_KEY'])
except KeyError as _:
    # No enviroment variable 'NCBI_API_KEY', try using '.netrc' file
    try:
        # Check .netrc exist and has no issues
        netrcFile = netrc.netrc()
        # Check there is an entry for NCBI api
        netrcFile.hosts['api.ncbi.nlm.nih.gov']
        # Use .netrc file by setting auth to none
    except netrc.NetrcParseError as e:
        print(e)
        print('WARNING: Please fix .netrc file before continuing')
        exit(-2)
    except FileNotFoundError as _:
        print("WARNING: No '.netrc' file and 'NCBI_API_KEY' environment variable was found.")
        print("WARNING: NCBI API request will be made WITHOUT authentication which will limit request rates")
    except KeyError as _:
        print("WARNING: No entry for 'api.ncbi.nlm.nih.gov' in  '.netrc' file and no 'NCBI_API_KEY' environment variable was found.")
        print("WARNING: NCBI API request will be made WITHOUT authentication which will limit request rates")


def _multi_attempt_get_request(url, max_attempts=5):
    '''
    This method will retry the get request if any issues occur.
    Arguments:
        url (str): the gest request url
        max_attempts (int, optional): the maximum number of retries allowed
    Returns:
        success (bool): If the request succeeded without exceeding max attempts
        resp : the get response or None if request exceeded max attempts

    '''
    attempts = 0
    while attempts < max_attempts:
        try:
            # Run request
            resp = get(url=url, auth=AUTH)
            # Check if we got a HTTP error
            resp.raise_for_status()
            # No errors return resp
            return True, resp
        except:
            attempts += 1
            print(f'Request failed, retrying (attempt {attempts}/{max_attempts})')
            sleep(1)
    print(f'Request failed, max attempts exceeded')
    return False, None

def get_genbank_dataset_reports_by_taxon(tax_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the genome dataset for dataset reports
    based on taxons. This request filters for full genomes, 
    that exactly match the given taxons from Genbank.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/taxon/-taxons-/dataset_report

    Arguments:
        tax_ids ([int]): A list of NCBI taxonmy ids
    
    Returns:
        A list of genome dataset reports
    '''

    if len(tax_ids) > 1000:
        raise RuntimeError('Too many taxonomy ids specified')
    url = (f'{GENOME_ENDPOINT}/taxon/{"%2C".join([str(x) for x in tax_ids])}' + 
            f'/dataset_report?filters.assembly_level=complete_genome&' + 
            f'filters.assembly_source=genbank&tax_exact_match=true&page_size={len(tax_ids)}')
    success, resp = _multi_attempt_get_request(url)
    return success, resp.json()['reports'] if success else resp

def get_refseq_dataset_reports_by_taxon(tax_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the genome dataset for dataset reports
    based on taxons. This request specifcally filters for full genomes, 
    that exactly match the given taxons from RefSeq.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/taxon/-taxons-/dataset_report

    Arguments:
        tax_ids ([int]): A list of NCBI taxonmy ids
    
    Returns:
        A list of genome dataset reports
    '''

    if len(tax_ids) > 1000:
        raise RuntimeError('Too many taxonomy ids specified')
    url = (f'{GENOME_ENDPOINT}/taxon/{"%2C".join([str(x) for x in tax_ids])}' + 
            f'/dataset_report?filters.assembly_level=complete_genome&' + 
            f'filters.assembly_source=refseq&tax_exact_match=true&page_size={len(tax_ids)}')
    success, resp = _multi_attempt_get_request(url)
    return success, resp.json()['reports'] if success else resp

def get_genes_by_id(gene_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, gets the gene fasta file for the given gene ids from the gene dataset.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/gene/id/-gene_ids-/download
    Arguments:
        gene_ids ([int]): A list of NCBI gene ids
    
    Returns:
        The zip file in bytes format
    '''
    url = (f'{GENE_ENDPOINT}/id/{"%2C".join([str(x) for x in gene_ids])}' +
            f'/download?include_annotation_type=FASTA_GENE')
    success, resp = _multi_attempt_get_request(url)
    return success, resp.content if success else resp

def get_assembly_by_accession(accessions):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, gets the fasta files for the given accessions from the genome dataset.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/accession/-accessions-/download
    Arguments:
        accessions ([int]): A list of NCBI assembly accession
    
    Returns:
        The zip file in bytes format
    '''
    url = (f'{GENOME_ENDPOINT}/accession/{"%2C".join(accessions)}' + 
            f'/download?include_annotation_type=GENOME_FASTA')
    success, resp = _multi_attempt_get_request(url)
    return success, resp.content if success else resp
