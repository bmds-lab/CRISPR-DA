'''
off-target.py

This file contains methods for assessing a CRISPR guide RNA's off-target risk against other genomes.
This is done using ISSL which is used in Crackling
Crackling can be found here: https://github.com/bmds-lab/Crackling
'''

import os
import multiprocessing
from pathlib import Path
from importlib import resources
from tempfile import TemporaryDirectory

from . import cache
from . import utils
from .collection import ViralCutCollection

def run_offtarget_scoring(collection: ViralCutCollection, accessions, processors=os.cpu_count()):
    '''Runs ISSL off-target scoring for the list of guides against each provided accession.

    Arguments:
        collection (ViralCut.Collection): A ViralCut collection containing the candidate guides.
        accessions (list): A list of accessions to be scored against.
        processors (int):  The number of processors to run. Zero indicates all available.

    Returns: None, all results are stored in the ViralCut Collection
    '''
    # Only process guides that were selected as efficient by on-target scoring
    guidesToScore = [g for g in collection if True in collection[g]['CDE_passed']]

    with resources.path('viralcut.resources', 'ISSLScoreOfftargets') as resource:
        isslScoreOfftargetsBin = resource

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
            args.append([f'{isslScoreOfftargetsBin} {isslIdx} {guidesFile} 4 0 and', resultsFile])

        errors = []
        # Begin scoring - single processor mode
        if processors == 1:
            for (command, stdOut), accession in zip(args,accessions):
                if not utils.run_command(command, stdOut):
                    errors.append(accession)
        # Begin scoring - multi-processor mode
        else:
            with multiprocessing.Pool(processors) as p:
                success = p.starmap(utils.run_command, args)
                errors = [accessions[idx] for idx, result in enumerate(success) if not result]

        # Add results to ViralCutCollection
        for accession in accessions:
            with open(resultsFiles[accession], 'r') as fp:
                lines = fp.readlines()
            lines = [line.strip().split('\t') for line in lines]
            for idx, [_, mit, cfd, uniqueSites, totalSites] in enumerate(lines):
                collection[guidesToScore[idx]].add_assembly_score(accession, mit, cfd, uniqueSites, totalSites)