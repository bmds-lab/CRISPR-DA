'''

Author: Jake Bradford

This script begins the project on scoring CRISPR-Cas9 guides against all viruses
in human.

Notes
    - NCBI Taxon ID for viruses is 10239
    - NCBI Taxon ID for Homo sapiens is 9606
    - NCBI Taxon ID for Homo sapiens neanderthalensis is 63221
    - NCBI Gene ID for SARS-CoV-2 spike or surface glycroprotein (S) gene is 43740568
    - NCBI Gene ID for SARS-CoV-2 envelope (E) gene is 43740570
'''

__all__ = ['run_analysis']

from . import data
from . import on_target
from . import off_target
from . import config
from . import visualise
from .collection import CRISPRDACollection

def run_analysis(target_accession = None, target_gene_id = None,  evaluation_accessions = None, evaluation_root_tax_id = None):
    '''Run pan-viral sgRNA design

    Arguments:
        target_accession (string):          The NCBI accession to target.
        target_gene_id (string):            The gene ID to extract sites from.
        evaluation_accessions (list):       A list of accessiion to evaluate
        evaluation_root_tax_id (string):    The taxon id for the root of the phylogentic tree to evaluate

    Returns:
        A CRISPRDACollection with results.
    '''
    if ((target_accession == None and target_gene_id == None) or 
        (target_accession != None and target_gene_id != None)) :
        print("Please provide either a target accession OR target gene id")
        exit(-1)

    if ((evaluation_accessions == None and evaluation_root_tax_id == None) or 
        (evaluation_accessions != None and evaluation_root_tax_id != None)) :
        print("Please provide either a list of accessions OR a taxon id for the root of the phylogentic tree to evaluate against")
        exit(-1)

    if config.VERBOSE:
        print('Downloading target sequence')
    if target_gene_id:
        target = data.download_ncbi_genes([target_gene_id])
    elif target_accession:
        target = data.download_ncbi_assemblies([target_accession])
    # Returns list (should get if one or more files, should one)
    if len(target) > 1:
        raise RuntimeError('There should only be one target sequence')
    target = target[0]

    if config.VERBOSE:
        print('Extracting guides')
    guides = data.extract_guides(target)

    if config.VERBOSE:
        print('Creating CRISPRDACollection')
    collection = data.create_collection(target, guides)

    if evaluation_root_tax_id:
        evaluation_tax_ids, evaluation_accessions = data.get_accession_from_root(evaluation_root_tax_id)

    collection.accessions =  evaluation_accessions
    if config.VERBOSE:
        print('Downloading assemblies')
    data.download_ncbi_assemblies(collection.accessions)

    if evaluation_accessions:
        evaluation_tax_ids = data.get_tax_ids_from_accessions(evaluation_accessions)

    collection.update_phylogenetic_tree(evaluation_tax_ids)
    collection.accession_to_tax_id = {
        accs : tax_id
        for accs, tax_id in zip(evaluation_accessions, evaluation_tax_ids)
    }

    # Create ISSL indexes for the downloaded accessions
    if config.VERBOSE:
        print('Generating ISSL indexes')
    data.create_issl_indexes(evaluation_accessions)

    # Evaluate on-target efficiency via Crackling    
    if config.VERBOSE:
        print('Evaluating efficiency')
    on_target.run_CRISPR_DeepEnsemble(collection, uncertainty_threshold=0.90)
    on_target.run_crackling_on_target(collection)

    # Do off-target scoring
    if config.VERBOSE:
        print('Assessing off-target risk')
    off_target.run_offtarget_scoring(collection, evaluation_accessions)

    # Calculate node scores
    collection.calculate_node_scores()

    visualise.generate_itol_tree(collection)

    if config.VERBOSE:
        print('Done.')
    
    return collection

