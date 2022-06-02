import argparse

import pandas as pd

from viralcut.score import run_analysis

def main():
    parser = argparse.ArgumentParser(
        description='ViralCut: for designing CRISPR-Cas9 sgRNA when many viral genomes are to be considered.'
    )
    
    parser.add_argument('-g', '--geneid', required=True, 
        help='The NCBI gene ID to extract potential CRISPR target sites from')
    parser.add_argument('-a', '--accessions', required=True, nargs='+', 
    help='A list of NCBI accessions to score CRISPR sites against')
    parser.add_argument('-o', '--output', required=False, default=None,
        help='Output files prefix (e.g. `/home/user/results/covid` would become `/home/user/results/covid-guides.csv`, etc.). Default is `--geneid`.')
    
    args = parser.parse_args()
    
    if args.output is None:
        args.output = args.geneid
    
    gc, scores = run_analysis(args.geneid, args.accessions)
    
    with open(f"{args.output}-guides.csv", "w") as fp:
        gc.to_dataframe().to_csv(fp, index=False)
            
    with open(f"{args.output}-scores.csv", "w") as fp:
        dfScores = pd.DataFrame(scores)
        dfScores['accession'] = dfScores.index
        accs_col = dfScores.pop('accession')
        dfScores.insert(0, 'accession', accs_col)
        dfScores.to_csv(fp, index=False)

    with open(f"{args.output}-guides.md", "w") as fp:
        gc.to_dataframe().to_markdown(fp, index=False)
            
    with open(f"{args.output}-scores.md", "w") as fp:
        dfScores = pd.DataFrame(scores)
        dfScores['accession'] = dfScores.index
        accs_col = dfScores.pop('accession')
        dfScores.insert(0, 'accession', accs_col)
        dfScores.to_markdown(fp, index=False)
            
    

if __name__ == '__main__':
    main()
