'''

Author: Jake Bradford

python3.8 design_for_covid_spike_protein.py

Notes
    - NCBI Gene ID for SARS-CoV-2 spike or surface glycroprotein (S) gene is 43740568
'''

from viralcut import run_analysis

import pandas as pd

def main():
    number_viruses = 50
    gene_id = '43740568'

    viruses = pd.read_csv('viruses.csv')
    accessions = list(viruses['Assembly'].head(number_viruses))

    gc, scores = run_analysis(target_gene_id=gene_id, evaluation_accessions=accessions)
    
    print(gc.to_dataframe())
    print(pd.DataFrame(scores))

if __name__ == '__main__':
    main()