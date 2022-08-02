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

from collections import namedtuple
from glob import glob
from tempfile import NamedTemporaryFile


from . import config
from . import analysis
from .design import process_gene_by_id, run_mini_crackling
from .data import *


def _run_issl_score(accession, path_guides, path_stdout):
    # assume the index exists and there is only one.
    issl_idx = glob(os.path.join(get_assembly_cache(accession), "*.issl"))[0]
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

    # Key: accession, Value: temporary file path
    issl_output_files = {}

    # Begin scoring - single processor mode
    if processors == 1:

        for accs in accessions:
            fpStdout = NamedTemporaryFile(delete=False)
            _run_issl_score(accession, guides_file.name, fpStdout)
            issl_output_files[accs] = fpStdout.name

    # Begin scoring - multi-processor mode
    else:
        with multiprocessing.Pool(os.cpu_count() if not processors else processors) as p:

            args = []
            for accs in accessions:
                tmp_file = NamedTemporaryFile(delete=False).name
                issl_output_files[accs] = tmp_file
                args.append((accs, guides_file.name, tmp_file))

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
    for accs in issl_output_files:
        with open(issl_output_files[accs], 'r') as fp:
            for line in fp:
                seq, mit, cfd, uniqueSites, totalSites, *_ = [x.strip() for x in line.split('\t')]

                scores['accession'].append(accs)
                scores['sequence'].append(seq)
                
                # reverse the global off-target score to what is now the 'assembly score'
                scores['mit'].append(10_000.0 / float(mit) - 100.0) 
                scores['cfd'].append(10_000.0 / float(cfd) - 100.0)
                
                scores['unique_sites'].append(uniqueSites)
                scores['total_sites'].append(totalSites)

    return scores


def run_analysis(gene_id, accessions=None):
    '''Run pan-viral sgRNA design

    Arguments:
        gene_id (string):   The gene ID to extract sites from.
        accessions (list):  (optional) A list of accessions to evaluate. If None then the 
                            pre-compiled list of human viruses will be used.

    Returns:
        A ViralCutCollection with results.
    '''
    
    if accessions is None:
        accessions = get_accessions_from_ncbi_table_export()

    # Download the requested accessions
    if config.VERBOSE:
        print('Downloading assemblies')

    accs_downloaded = download_ncbi_assemblies(accessions)

    # Create ISSL indexes for the downloaded accessions
    if config.VERBOSE:
        print('Generating ISSL indexes')

    create_issl_indexes(accessions)

    # Download the gene then extract CRISPR sites and score each
    if config.VERBOSE:
        print('Downloading gene sequence then analysing')

    download_ncbi_genes([gene_id])

    # Extract sites    
    if config.VERBOSE:
        print('Extracting target sites')
    collection = process_gene_by_id(gene_id)

    collection.gene_id = gene_id
    collection.accessions = accessions
    
    collection.accession_to_tax_id = {
        accs : tax_id
        for accs, tax_id in zip(accessions, analysis.get_tax_ids_from_accessions(accessions, uniq=False))
    }
    
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
    scores = run_offtarget_scoring(targets_to_score, accessions)

    # Add scores to the collection
    for score_name in config.ISSL_SCORES:
        for accession, seq, score, uniq, total in zip(
            scores['accession'],
            scores['sequence'],
            scores[score_name],
            scores['unique_sites'],
            scores['total_sites'],
        ):
            collection[seq].add_assembly_score(
                accession,
                score_name,
                score,
                uniq,
                total
            )

    # Calculate node scores
    #collection.calculate_node_scores()

    if config.VERBOSE:
        print('Done.')
    
    return collection

