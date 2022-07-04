from . import analysis
from .data import get_accessions_from_ncbi_table_export

from flask import Flask, send_from_directory, make_response

app = Flask(__name__,
    static_url_path='', 
    static_folder='viewer',
    template_folder='templates'
)

accessions = get_accessions_from_ncbi_table_export('/home/jake/ViralCut/examples/viruses.csv')[0:50]

tax_ids = analysis.get_tax_ids_from_accessions(accessions) + analysis.get_descendant_tax_ids_from_root_tax_id(root_tax_id=687331)

@app.route("/")
def root():
    return send_from_directory('viewer', 'index.html')

@app.route("/newick")
def newick():
    output = make_response(
        analysis.generate_newick_string_from_tax_ids(tax_ids)
    )
    output.headers["Content-type"] = "text/plain"
    return output
    

@app.route("/scores")
def scores():
    df = analysis.generate_df_phylo_node_scores_from_tax_ids(tax_ids)
    df['score'] = df['score'].round(decimals=2)
    output = make_response(
        df.to_csv(index=False)
    )
    output.headers["Content-type"] = "text/csv"
    return output

@app.route("/species")
def species():
    print(tax_ids)
    df = analysis.generate_df_of_species_data_from_tax_ids(tax_ids)
    output = make_response(
        df.to_csv(index=False)
    )
    output.headers["Content-type"] = "text/csv"
    return output

def run():
    app.run(host='0.0.0.0', port='8080', debug=True)

if __name__ == '__main__':
    run()
    
    
    
    
    
    
    

def generate_df_phylo_node_scores_from_accessions(accessions):
    '''Given a list of accessions and scores, calculate node scores.
    
    Arguments:
        accessions (list): A list of accessions
        
    Returns:
        A DataFrame with columns: tax_id, score
    
    '''
    with open('/home/jake/ViralCut/43740568-scores.csv', 'r') as fp:
        df = pd.read_csv(fp)
    
    df['local'] = 10000.0 / df['mit'] - 100.0
    
    tax_ids = get_tax_ids_from_accessions(accessions)
    print(tax_ids)
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(tax_ids)#, intermediate_nodes=True)

    data = {'tax_id': [], 'score' : []}

    for idx, i in enumerate(tree.traverse(strategy="levelorder")):
        species = map(float, i.get_leaf_names())
        df_species = df[df['taxId'].isin(species)]
        mit = 10000.0 / (100.0 + df_species['local'].sum())
        
        data['tax_id'].append(i.name)
        data['score'].append(mit)
    
    return pd.DataFrame(data)