import os
import netrc
from requests import get
from requests.auth import HTTPBasicAuth

GENOME_ENDPOINT = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome'
GENE_ENDPOINT = 'https://api.ncbi.nlm.nih.gov/datasets/v2gene'
AUTH = None

# Check for NCBI API key
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


def get_genbank_dataset_reports_by_taxon(tax_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the genome dataset for dataset reports
    based on taxons. This request specifcally filters for full 
    genomes, that exactly match the given taxons from Genbank.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/taxon/-taxons-/dataset_report

    Arguments:
        tax_ids ([int]): A list of NCBI taxonmy ids
    
    Returns:
        A list of genome dataset reports
    '''

    if len(tax_ids) > 1000:
        raise RuntimeError('Too many taxonomy ids specified')
    # Query NCBI dataset API
    resp = get(url=f'{GENOME_ENDPOINT}/taxon/{"%2C".join([str(x) for x in tax_ids])}/dataset_report?filters.assembly_level=complete_genome&filters.assembly_source=genbank&tax_exact_match=true&page_size={len(tax_ids)}', auth=AUTH)
    # Check if we got a HTTP error
    resp.raise_for_status()
    # Response should be in json format
    resp = resp.json()
    return resp['reports']

def get_refseq_dataset_reports_by_taxon(tax_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the genome dataset for dataset reports
    based on taxons. This request specifcally filters for full 
    genomes, that exactly match the given taxons from RefSeq.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/taxon/-taxons-/dataset_report

    Arguments:
        tax_ids ([int]): A list of NCBI taxonmy ids
    
    Returns:
        A list of genome dataset reports
    '''

    if len(tax_ids) > 1000:
        raise RuntimeError('Too many taxonomy ids specified')
    resp = get(url=f'{GENOME_ENDPOINT}/taxon/{"%2C".join([str(x) for x in tax_ids])}/dataset_report?filters.assembly_level=complete_genome&filters.assembly_source=refseq&tax_exact_match=true&page_size={len(tax_ids)}', auth=AUTH)
    # Check if we got a HTTP error
    resp.raise_for_status()
    resp = resp.json()
    return resp['reports']

def get_genes_by_id(gene_ids):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the gene dataset for the gene fasta file.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/gene/id/-gene_ids-/download
    Arguments:
        gene_ids ([int]): A list of NCBI gene ids
    
    Returns:
        The zip file in bytes format
    '''
    # Query NCBI dataset API
    resp = get(url=f'{GENE_ENDPOINT}/id/{"%2C".join(gene_ids)}/download?include_annotation_type=FASTA_GENE', auth=AUTH)
    # Check if we got a HTTP error
    resp.raise_for_status()
    return resp.content

def get_assembly_by_accession(accessions):
    '''
    This functions queries the NCBI Dataset V2 REST API.
    Specifcally, queries the genome dataset for the accessions fasta file.
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#get-/genome/accession/-accessions-/download
    Arguments:
        accessions ([int]): A list of NCBI assembly accession
    
    Returns:
        The zip file in bytes format
    '''
    # Query NCBI dataset API
    resp = get(url=f'{GENOME_ENDPOINT}/accession/{"%2C".join(accessions)}/download?include_annotation_type=GENOME_FASTA', auth=AUTH)
    # Check if we got a HTTP error
    resp.raise_for_status()
    return resp.content
