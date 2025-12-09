'''
on-target.py

This file contains methods for assessing a CRISPR guide RNA's on-target effiency.
CRISPR DeepEnsemble can be found here: https://github.com/bmds-lab/CRISPR_DeepEnsemble 
Crackling can be found here: https://github.com/bmds-lab/Crackling/blob/9e9d78196e97fe11e60f0d9bcc7c7e1349a03ae4/src/crackling/Crackling.py
'''

import sys
import ast
import os
import re
import subprocess
import joblib
import importlib
import shutil
import numpy as np
import torch as t
import pandas as pd
from importlib import resources 
from tempfile import NamedTemporaryFile
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from . import utils
from . import CRISPR_DeepEnsemble
sys.modules['CRISPR_DeepEnsemble'] = CRISPR_DeepEnsemble

CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNKNOWN = '?'
CODE_ERROR = '!'

config = {}
config['rnafold'] = {}
config['sgrnascorer2'] = {}
config['rnafold']['output'] = NamedTemporaryFile(delete=False).name
config['rnafold']['input'] = NamedTemporaryFile(delete=False).name
config['rnafold']['binary'] = 'RNAfold'
config['rnafold']['threads'] = 32
config['rnafold']['low_energy_threshold'] = -30
config['rnafold']['high_energy_threshold'] = -18
config['sgrnascorer2']['model'] = 'model-1_0_2.txt'
config['sgrnascorer2']['score_threshold'] = 0

t.manual_seed(123)
dtype = t.float64
t.set_default_dtype(dtype)

def run_CRISPR_DeepEnsemble(collection, score_threshold=0.7, uncertainty_threshold=0.05):
    '''
    Utilising the model from https://github.com/bmds-lab/CRISPR_DeepEnsemble to predict on-target efficiency.
    'This allows to design strategies that consider uncertainty in guide RNA selection'.

    Arguments:
    collection (CRISPRDA.Collection): A CRISPRDA collection containing the candidate guides.
    score_threshold (float): A score threshold that must be exceed for a candidate guide to be 'efficient'. 
                             Must be in the range from 0 to 1.
    uncertainty_threshold (float): A uncertainty threshold that must not be exceed for a score to be considered 'acceptable'. 
                                   Must be in the range from 0 to 1.

    Returns: None, all results are stored in the CRISPRDA Collection
    '''
    # NOTE
    # CDE = CRISPR Deep Ensemble
    # UQ = Uncertainty Quantification 

    # Load deep ensemble from resources
    with resources.path('viralcut.resources', 'CRISPR_DeepEnsemble.zip') as model:
        ensemble = CRISPR_DeepEnsemble.RegressionDeepEnsemble(load_from=model)

    # Convert threshold percent to threshold value
    with resources.path('viralcut.resources', 'trainingResults.pkl') as trainingResults:
        trainingResults = pd.read_pickle(trainingResults)
    UQ_threshold = np.quantile(trainingResults["range"].to_numpy(), [uncertainty_threshold], interpolation="nearest")

    # Encode data and extract features
    oneHot = []
    for guide in collection:
        for target30 in collection[guide]['30mer']:
            oneHot.append(np.array(utils.one_hot_encode(target30)).tolist())
    meltingPoint = []
    for guide in collection:
        for target30 in collection[guide]['30mer']:
            myseq = Seq(target30)
            meltingPoint.append(mt.Tm_NN(myseq))
    _onehot = t.tensor(oneHot).transpose(1,2).unsqueeze(dim=1) 
    _meltingpoint = t.tensor(meltingPoint).reshape(-1,1)

    # Score guides
    prediction = ensemble.predict(inputs = (_onehot, _meltingpoint)).tolist()
    UQLowerBound, UQUpperBound, UQInterquartileRange = [x.tolist() for x in ensemble.uncertainty_bounds(inputs = (_onehot, _meltingpoint), n_samples=1000, lower=0.01, upper=0.99)]

    # Add scores to CRISPRDA Collection
    for guide in collection:
        for target30 in collection[guide]['30mer']:
            score = prediction.pop(0)
            UQIQR = UQInterquartileRange.pop(0)
            UQRange = UQUpperBound.pop(0) - UQLowerBound.pop(0)
            if hasattr(collection[guide], 'CDE_score'):
                collection[guide]['CDE_score'].append(score)
                collection[guide]['CDE_UQ_range'].append(UQRange)
                collection[guide]['CDE_UQ_IQR'].append(UQIQR)
                if score > score_threshold and UQRange < UQ_threshold:
                    collection[guide]['CDE_passed'].append(True)
                else:
                    collection[guide]['CDE_passed'].append(False)
            else:
                collection[guide]['CDE_score'] = [score]
                collection[guide]['CDE_UQ_range'] = [UQRange]
                collection[guide]['CDE_UQ_IQR'] = [UQIQR]
                if score > score_threshold and UQRange < UQ_threshold:
                    collection[guide]['CDE_passed'] = [True]
                else:
                    collection[guide]['CDE_passed'] = [False]



def run_crackling_on_target(collection):
    '''
    A subset of the Crackling method that scores the on-target efficiency.
    Based on: https://github.com/bmds-lab/Crackling/blob/9e9d78196e97fe11e60f0d9bcc7c7e1349a03ae4/src/crackling/Crackling.py.

    Arguments:
        collection (CRISPRDA.Collection): A CRISPRDA collection containing the candidate guides.
    
    Returns: None, all results are stored in the CRISPRDA Collection
    '''

    ## G20
    for target23 in collection:
        if target23[19] != 'G':
            collection[target23]['passed_g20'] = CODE_REJECTED
        else:
            collection[target23]['passed_g20'] = CODE_ACCEPTED

    ## Removing targets with leading T
    for target23 in collection:
        if (target23[-2:] == 'GG' and target23[0] == 'T') or \
            (target23[:2] == 'CC' and target23[-1] == 'A'):
            collection[target23]['passed_avoid_leading_t'] = CODE_REJECTED
        else:
            collection[target23]['passed_avoid_leading_t'] = CODE_ACCEPTED

    ## AT% ideally is between 20-65%
    for target23 in collection:
        at_percentage = float(sum([x.upper() in ['A', 'T'] for x in target23[0:20]]) / 20.0 * 100.0)

        if at_percentage < 20 or at_percentage > 65:
            collection[target23]['passed_at_percent'] = CODE_REJECTED
        else:
            collection[target23]['passed_at_percent'] = CODE_ACCEPTED

        collection[target23]['AT'] = at_percentage

    ## Removing targets that contain TTTT
    for target23 in collection:
        if 'TTTT' in target23:
            collection[target23]['passed_tttt'] = CODE_REJECTED
        else:
            collection[target23]['passed_tttt'] = CODE_ACCEPTED

    ## Calculating secondary structures
    guide = 'GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU'
    pattern_rna_structure = r'.{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)'
    pattern_rna_energy = r'\s\((.+)\)'

    if os.path.exists(config['rnafold']['output']):
        os.remove(config['rnafold']['output'])

    with open(config['rnafold']['input'], 'w+') as f_rna_input:
        for target23 in collection:
            f_rna_input.write(f'G{target23[1:20]}{guide}\n')

    subprocess.run(
        map(str, [
            config['rnafold']['binary'],
            '--noPS',
            f"-j{config['rnafold']['threads']}",
            f"-i{config['rnafold']['input']}",
            '-o'
        ]),
        #shell=True,
        check=True
    )

    shutil.move('RNAfold_output.fold', config['rnafold']['output'])
    rna_structures = {}
    with open(config['rnafold']['output'], 'r') as fRnaOutput:
        i = 0
        line_1, line_2, target = None, None, None
        for line in fRnaOutput:
            if i % 2 == 0:
                # 0th, 2nd, 4th, etc.
                line_1 = line.rstrip()
                target = line_1[0:20]
            else:
                # 1st, 3rd, 5th, etc.
                line_2 = line.rstrip()
                rna_structures[utils.trans_to_dna(target[1:20])] = [
                    line_1, line_2, target
                ]

            i += 1


    for target23 in collection:
        key = target23[1:20]

        if key not in rna_structures:
            continue

        line_1 = rna_structures[key][0]
        line_2 = rna_structures[key][1]
        target = rna_structures[key][2]

        structure = line_2.split(' ')[0]
        energy = line_2.split(' ')[1][1:-1]

        collection[target23]['ss_line_1'] = line_1
        collection[target23]['ss_structure'] = structure
        collection[target23]['ss_energy'] = energy

        if utils.trans_to_dna(target) != target23[0:20] and utils.trans_to_dna('C'+target[1:]) != target23[0:20] and utils.trans_to_dna('A'+target[1:]) != target23[0:20]:
            collection[target23]['passed_secondary_structure'] = CODE_ERROR
            continue

        match_structure = re.search(pattern_rna_structure, line_2)
        if match_structure:
            energy = ast.literal_eval(match_structure.group(1))
            if energy < float(config['rnafold']['low_energy_threshold']):
                collection[utils.trans_to_dna(target23)]['passed_secondary_structure'] = CODE_REJECTED
            else:
                collection[target23]['passed_secondary_structure'] = CODE_ACCEPTED
        else:
            match_energy = re.search(pattern_rna_energy, line_2)
            if match_energy:
                energy = ast.literal_eval(match_energy.group(1))
                if energy <= float(config['rnafold']['high_energy_threshold']):
                    collection[utils.trans_to_dna(target23)]['passed_secondary_structure'] = CODE_REJECTED
                else:
                    collection[target23]['passed_secondary_structure'] = CODE_ACCEPTED

    ## Calc mm10db result
    for target23 in collection:
        if not all([
            collection[target23]['passed_at_percent'] == CODE_ACCEPTED,
            collection[target23]['passed_tttt'] == CODE_ACCEPTED,
            collection[target23]['passed_secondary_structure'] == CODE_ACCEPTED,
            collection[target23]['passed_avoid_leading_t'] == CODE_ACCEPTED,
        ]):
            collection[target23]['accepted_by_mm10db'] = CODE_REJECTED
        else:
            collection[target23]['accepted_by_mm10db'] = CODE_ACCEPTED

    ## sgRNAScorer 2.0 model
    encoding = {
        'A' : '0001',    'C' : '0010',    'T' : '0100',    'G' : '1000',
        'K' : '1100',    'M' : '0011',    'R' : '1001',    'Y' : '0110',
        'S' : '1010',    'W' : '0101',    'B' : '1110',    'V' : '1011',
        'H' : '0111',    'D' : '1101',    'N' : '1111'
    }

    with importlib.resources.path('viralcut.resources', config['sgrnascorer2']['model']) as fp:
        clf_linear = joblib.load(fp) #config['sgrnascorer2']['model'])

    for target23 in collection:
        sequence = target23.upper()
        entry_list = []

        for x in range(0, 20):
            for y in range(0, 4):
                entry_list.append(int(encoding[sequence[x]][y]))

        # predict based on the entry
        score = clf_linear.decision_function([entry_list])[0]

        collection[target23]['sgrnascorer2_score'] = score

        if float(score) < float(float(config['sgrnascorer2']['score_threshold'])):
            collection[target23]['accepted_by_sgrnascorer'] = CODE_REJECTED
        else:
            collection[target23]['accepted_by_sgrnascorer'] = CODE_ACCEPTED

    ## Begin efficacy consensus
    for target23 in collection:
        collection[target23]['consensus_count'] = sum([
            collection[target23]['accepted_by_mm10db'] == CODE_ACCEPTED,
            collection[target23]['accepted_by_sgrnascorer'] == CODE_ACCEPTED,
            collection[target23]['passed_g20'] == CODE_ACCEPTED,
        ])
        collection[target23]['passed_crackling'] = CODE_ACCEPTED if collection[target23]['consensus_count'] >= 2 else CODE_REJECTED
