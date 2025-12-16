'''

Author: Jake Bradford

python3.8 design_for_covid_spike_protein.py
A minimal example of how to access CRISPR-DA programmatically

Notes
    - NCBI Gene ID for SARS-CoV-2 spike or surface glycroprotein (S) gene is 43740568
'''

from crispr_da import run_analysis
import pandas as pd

def main():
    number_viruses = 50
    gene_id = '43740568'

    # NOTE: This was used for testing purposes, Please update later
    viruses = pd.read_csv('viruses.csv')
    accessions = list(viruses['Assembly Accession'].head(number_viruses))

    run_analysis(target_gene_id=gene_id, evaluation_accessions=accessions)

if __name__ == '__main__':
    main()