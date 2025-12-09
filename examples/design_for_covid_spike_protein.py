'''

Author: Jake Bradford

python3.8 design_for_covid_spike_protein.py

Notes
    - NCBI Gene ID for SARS-CoV-2 spike or surface glycroprotein (S) gene is 43740568
'''

from crisprda import run_analysis

import pandas as pd

def main():
    number_viruses = 50
    gene_id = '43740568'

    # NOTE: This was used for testing purposes, Please update later
    viruses = pd.read_csv('viruses_updated.csv')
    accessions = list(viruses['Assembly Accession'].head(number_viruses))
    output = run_analysis(target_gene_id=gene_id, evaluation_accessions=accessions)
    g = []
    a  = []
    mit = []
    cfd = []
    total = []
    unique = []
    for accession in accessions:
        for guide in output.guides:
            if len(output.guides[guide].assembly_scores) > 0:
                m, c, unique_sites, total_sites = output.guides[guide].assembly_scores[accession].values()
                g.append(guide)
                a.append(accession)
                mit.append(m)
                cfd.append(c)
                total.append(total_sites)
                unique.append(unique_sites)

    data = pd.DataFrame({'guide': g, 'accession': a, 'MIT score': mit, 'CFD score': cfd, 'total_sites': total, 'unique_sites': unique})
    # data.to_csv('/mnt/ssd1/carl/ViralCut/issl_scores.csv', index=False)

    # print(gc.to_dataframe())
    # print(pd.DataFrame(scores))

if __name__ == '__main__':
    main()