'''
off-target.py

This file contains methods for assessing a CRISPR guide RNA's off-target risk against other genomes.
This is done using ISSL which is used in Crackling
Crackling can be found here: https://github.com/bmds-lab/Crackling
'''

import os
import subprocess
import multiprocessing
from pathlib import Path
from tempfile import TemporaryDirectory

from . import cache
from . import config
from .collection import ViralCutCollection

def _run_issl_score(path_guides, path_issl_index, path_results):
    '''
    Runs ISSL off-target scoring with the given arguments.
    This is used to run scoring using mulitple processors.

    Arguments:
        path_guides (path_like str): The file that contains the guides to be scored
        path_issl_index (path_like str): The ISSL index file
        path_results (path_like str): The output file, where the results will be stored
    
    Returns: None
    '''
    # assume the index exists and there is only one.
    with open(path_results, 'w') as outFile:
        subprocess.run([
            config.BIN_ISSL_SCORE,
            path_issl_index,
            path_guides,
            '4',
            '0',
            'and'
        ], stdout=outFile)


def run_offtarget_scoring(collection: ViralCutCollection, accessions, processors=0):
    '''Runs ISSL off-target scoring for the list of guides against each provided accession.

    Arguments:
        collection (ViralCut.Collection): A ViralCut collection containing the candidate guides.
        accessions (list): A list of accessions to be scored against.
        processors (int):  The number of processors to run. Zero indicates all available.

    Returns: None, all results are stored in the ViralCut Collection
    '''
    # Only process guides that were selected as efficient by on-target scoring
    guidesToScore = [g for g in collection if True in collection[g]['CDE_passed']]

    # Create temporary working dir
    with TemporaryDirectory() as tmpDir:
        tmpPath = Path(tmpDir)
        # Write the guides that passed on-target scoring to a temporary file
        guidesFile = tmpPath / 'guides.txt'
        with open(guidesFile, 'w') as outFile:
            for guide in guidesToScore:
                outFile.write(f"{guide[:20]}\n")

        # Key: accession, Value: temporary file path
        resultsFiles = {}
        # For running scoring in multiprocessing mode
        args = []
        for accession in accessions:
            isslIdx = cache.get_file(accession, '.issl')
            resultsFile = tmpPath / f'results_{accession}.txt'
            resultsFiles[accession] = resultsFile
            args.append((guidesFile, isslIdx, resultsFile))

        # Begin scoring - single processor mode
        if processors == 1:
            for guideFile, isslIdx,  resultsFile in args:
                _run_issl_score(guideFile, isslIdx, resultsFile)
        # Begin scoring - multi-processor mode
        else:
            with multiprocessing.Pool(os.cpu_count() if not processors else processors) as p:
                p.starmap(_run_issl_score, args)

        # Add results to ViralCutCollection
        for accession in accessions:
            with open(resultsFiles[accession], 'r') as fp:
                lines = fp.readlines()
            lines = [line.strip().split('\t') for line in lines]
            for idx, [_, mit, cfd, uniqueSites, totalSites] in enumerate(lines):
                collection[guidesToScore[idx]].add_assembly_score(accession, mit, cfd, uniqueSites, totalSites)