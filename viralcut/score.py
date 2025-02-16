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

import os
import subprocess
import multiprocessing
from tempfile import NamedTemporaryFile


from . import config
from .design import run_mini_crackling, run_crispr_deep_ensemble
from .data import *

def _run_issl_score(accession, path_guides, path_stdout):
    # assume the index exists and there is only one.
    issl_idx = cache.get_accession_issl_index(accession)
    with open(path_stdout, 'w') as fp:
        subprocess.run([
            config.BIN_ISSL_SCORE,
            issl_idx,
            path_guides,
            '4',
            '0',
            'and'
        ], stdout=fp)


def run_offtarget_scoring(guides, accessions, processors=0):
    '''Runs ISSL off-target scoring for the list of guides against each provided accession.

    Arguments:
        guides (list):     A list of guides to score
        accessions (list): A list of accessions to be scored against.
        processors (int):  The number of processors to run. Zero indicates all available.

    Returns:
        Nested dictionaries:
            Level-1 key: guide
            Level-2 key: accession
            Level-2 val: vector of MIT then CFD scores
    '''

    # Write the guides to a temporary file
    guides_file = NamedTemporaryFile(delete=False)
    with open(guides_file.name, 'w') as fp:
        for guide in guides:
            fp.write(f"{guide[:20]}\n")

    # guides_file = '/mnt/ssd1/carl/ViralCut/examples/genomes-viruses/data/issl_guides.txt'
    # Key: accession, Value: temporary file path
    issl_output_files = {}

    # Begin scoring - single processor mode
    if processors == 1:
        for accs in accessions:
            fpStdout = NamedTemporaryFile(delete=False).name
            # _run_issl_score(accession, guides_file, fpStdout)
            _run_issl_score(accs, guides_file.name, fpStdout)
            issl_output_files[accs] = fpStdout

    # Begin scoring - multi-processor mode
    else:
        args = []
        for accs in accessions:
            tmp_file = NamedTemporaryFile(delete=False).name
            issl_output_files[accs] = tmp_file
            # args.append((accs, guides_file, tmp_file))
            args.append((accs, guides_file.name, tmp_file))
        with multiprocessing.Pool(os.cpu_count() if not processors else processors) as p:
            p.starmap(_run_issl_score, args)

    # Collect all the scores
    scores = {
        'accession' :   [],
        'sequence' :    [],
        'mit' :         [],
        'cfd' :         [],
        'unique_sites' : [],
        'total_sites' :  [],
    }

    # print(issl_output_files)
    for accs in issl_output_files:
        # shutil.move(issl_output_files[accs], f'/mnt/ssd1/carl/ViralCut/examples/genomes-viruses/data/issl_offtargets/{accs}.txt')
        with open(issl_output_files[accs], 'r') as fp:
            lines = fp.readlines()
        lines = [line.strip().split('\t') for line in lines]

        for idx, (_, mit, cfd, uniqueSites, totalSites) in enumerate(lines):
            scores['accession'].append(accs)
            scores['sequence'].append(guides[idx])
            
            # reverse the global off-target score to what is now the 'assembly score'
            # scores['mit'].append(10_000.0 / float(mit) - 100.0) 
            # scores['cfd'].append(10_000.0 / float(cfd) - 100.0)
            scores['mit'].append(float(mit)) 
            scores['cfd'].append(float(cfd))

            scores['unique_sites'].append(uniqueSites)
            scores['total_sites'].append(totalSites)

    return scores

def process_target_accession(target_accession):
    if config.VERBOSE:
        print('Downloading target sequence')
    download_ncbi_assemblies([target_accession])
    if config.VERBOSE:
        print('Extracting guides')
    guides = get_guides_from_genome(target_accession)
    if config.VERBOSE:
        print('Creating ViralCutCollection')
    collection = create_collection_from_accession(target_accession, guides)
    return collection

def process_target_gene(target_gene_id):
    if config.VERBOSE:
        print('Downloading target sequence')
    download_ncbi_genes([target_gene_id])
    if config.VERBOSE:
        print('Extracting guides')
    guides = get_guides_from_gene(target_gene_id)
    if config.VERBOSE:
        print('Creating ViralCutCollection')
    collection = create_collection_from_gene(target_gene_id, guides)
    return collection

def process_evaluation_accessions(collection: ViralCutCollection, evaluation_accessions):
    collection.accessions =  evaluation_accessions
    evaluation_tax_ids = get_tax_ids_from_accessions(evaluation_accessions)
    collection.update_phylogenetic_tree(evaluation_tax_ids)
    collection.accession_to_tax_id = {
        accs : tax_id
        for accs, tax_id in zip(evaluation_accessions, evaluation_tax_ids)
    }
    return evaluation_accessions

def process_evaluation_tax_id(collection: ViralCutCollection, evaluation_root_tax_id):
    collection.reset_phylogenetic_tree(evaluation_root_tax_id)
    evaluation_tax_ids = [node['taxid'] for node in collection._ncbi_tree.traverse("postorder")]
    evaluation_accessions = get_accession_from_tax_id(evaluation_tax_ids)
    collection.accessions =  evaluation_accessions
    collection.accession_to_tax_id = {
        accs : tax_id
        for accs, tax_id in zip(evaluation_accessions, evaluation_tax_ids)
    }
    return evaluation_accessions


def run_analysis(target_accession = None, target_gene_id = None,  evaluation_accessions = None, evaluation_root_tax_id = None):
    '''Run pan-viral sgRNA design

    Arguments:
        target_accession (string):          The NCBI accession to target.
        target_gene_id (string):            The gene ID to extract sites from.
        evaluation_accessions (list):       A list of accessiion to evaluate
        evaluation_root_tax_id (string):    The taxon id for the root of the phylogentic tree to evaluate

    Returns:
        A ViralCutCollection with results.
    '''
    if ((target_accession == None and target_gene_id == None) or 
        (target_accession != None and target_gene_id != None)) :
        print("Please provide either a target accession OR target gene id")
        exit(-1)

    if ((evaluation_accessions == None and evaluation_root_tax_id == None) or 
        (evaluation_accessions != None and evaluation_root_tax_id != None)) :
        print("Please provide either a list of accessions OR a taxon id for the root of the phylogentic tree to evaluate against")
        exit(-1)

    if target_gene_id:
        collection = process_target_gene(target_gene_id)
    elif target_accession:
        collection = process_target_accession(target_accession)

    if evaluation_accessions:
        evaluation_accessions = process_evaluation_accessions(collection, evaluation_accessions)
    elif evaluation_root_tax_id:
        evaluation_accessions = process_evaluation_tax_id(collection, evaluation_root_tax_id)

    # Download the requested accessions
    if config.VERBOSE:
        print('Downloading assemblies')
    download_ncbi_assemblies(evaluation_accessions)

    # Create ISSL indexes for the downloaded accessions
    if config.VERBOSE:
        print('Generating ISSL indexes')
    create_issl_indexes(evaluation_accessions)

    # Evaluate on-target efficiency via Crackling    
    if config.VERBOSE:
        print('Evaluating efficiency')
    run_mini_crackling(collection)

    # Filter out guides that should not be assessed for off-target risk
    if config.VERBOSE:
        print('Assessing off-target risk')
        
    targets_to_score = []
    for guide in collection:
        if collection[guide]['consensus_count'] >= config.CONSENSUS_N:
            targets_to_score.append(guide)

    # Do off-target scoring
    scores = run_offtarget_scoring(targets_to_score, evaluation_accessions)

    # Add scores to the collection
    # for score_name in config.ISSL_SCORES:
    for accession, seq, mit, cfd, uniq, total in zip(
        scores['accession'],
        scores['sequence'],
        scores['mit'],
        scores['cfd'],
        scores['unique_sites'],
        scores['total_sites'],
    ):
        collection[seq].add_assembly_score(
            accession,
            mit,
            cfd,
            uniq,
            total
        )

    # Calculate node scores
    collection.calculate_node_scores()

    if config.VERBOSE:
        print('Done.')
    
    return collection

