'''
on-target.py

This file contains methods for assessing a CRISPR guide RNA's on-target effiency.
CRISPR DeepEnsemble can be found here: https://github.com/bmds-lab/CRISPR_DeepEnsemble 
'''

import sys
import numpy as np
import torch as t
import pandas as pd
from importlib import resources 
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from . import utils
from . import CRISPR_DeepEnsemble
sys.modules['CRISPR_DeepEnsemble'] = CRISPR_DeepEnsemble

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
    with resources.path('crisprda.resources', 'CRISPR_DeepEnsemble.zip') as model:
        ensemble = CRISPR_DeepEnsemble.RegressionDeepEnsemble(load_from=model)

    # Convert threshold percent to threshold value
    with resources.path('crisprda.resources', 'trainingResults.pkl') as trainingResults:
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

