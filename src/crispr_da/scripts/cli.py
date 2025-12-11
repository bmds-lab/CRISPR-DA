import argparse
import importlib
from crispr_da import run_analysis
from crispr_da.config import run_config

def main():
    parser = argparse.ArgumentParser(prog="crispr_da", description='CRISPR-DA: for designing CRISPR-Cas9 sgRNA when many genomes are to be considered.')
    parser.add_argument('-v', '--version', help="Print CRISPR-DA version", action='version', version=f'%(prog)s version {importlib.metadata.version('crispr_da')}')
    subParsers = parser.add_subparsers(dest='command', title='subcommands')

    configParser = subParsers.add_parser('config', help='Run config')
    configParser.add_argument('-f', '--force', help='Force rebuild bins', action='store_true', dest='force')

    analysisParser = subParsers.add_parser('analyse', help='Run analysis')
    onTarget = analysisParser.add_mutually_exclusive_group(required=True)
    onTarget.add_argument('--target_accession', action='store', default=None, dest='target_accession')
    onTarget.add_argument('--target_gene_id', action='store', default=None, dest='target_gene_id')
    offTarget = analysisParser.add_mutually_exclusive_group(required=True)
    offTarget.add_argument('--evaluation_accessions', nargs='+', action='store', default=None, dest='evaluation_accessions')
    offTarget.add_argument('--evaluation_root_tax_id', action='store', default=None, dest='evaluation_root_tax_id')

    args = parser.parse_args()
    if args.command == 'config':
        run_config(args.force)
    elif args.command == 'analyse':
        run_analysis(args.target_accession, args.target_gene_id, args.evaluation_accessions, args.evaluation_root_tax_id)
    else:
        parser.print_help()
    

if __name__ == '__main__':
    main()
