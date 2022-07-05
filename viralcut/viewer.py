from . import analysis
from . import data

from flask import Flask, send_from_directory, make_response

app = Flask(__name__,
    static_url_path='', 
    static_folder='viewer',
    template_folder='templates'
)

MAX_TAX_IDS = 30

#accessions = data.get_accessions_from_ncbi_table_export('/home/jake/ViralCut/examples/viruses.csv')[0:50]
#tax_ids = analysis.get_tax_ids_from_accessions(accessions) + 

def get_tax_ids(root_tax_id):
    #return data.get_cached_tax_ids()[:MAX_TAX_IDS]
    tax_ids = analysis.get_descendant_tax_ids_from_root_tax_id(root_tax_id=root_tax_id, max_depth=3)
    return tax_ids[:MAX_TAX_IDS]
    

@app.route("/")
def root():
    return send_from_directory('viewer', 'index.html')

@app.route("/newick/<int:root_tax_id>")
def newick(root_tax_id):
    tax_ids = get_tax_ids(root_tax_id)
    output = make_response(
        analysis.generate_newick_string_from_tax_ids(tax_ids)
    )
    output.headers["Content-type"] = "text/plain"
    return output
    

@app.route("/scores/<int:root_tax_id>")
def scores(root_tax_id):
    tax_ids = get_tax_ids(root_tax_id)
    df = analysis.generate_df_phylo_node_scores_from_tax_ids(tax_ids)
    #df['score_orig'] = df['score'].round(decimals=2)
    #df['score'] = (df['score_orig'] - df['score_orig'].min()) / (df['score_orig'].max() - df['score_orig'].min()) * 100.0

    output = make_response(
        df.to_csv(index=False)
    )
    output.headers["Content-type"] = "text/csv"
    return output

@app.route("/species/<int:root_tax_id>")
def species(root_tax_id):
    tax_ids = get_tax_ids(root_tax_id)
    df = analysis.generate_df_of_species_data_from_tax_ids(tax_ids)
    output = make_response(
        df.to_csv(index=False)
    )
    output.headers["Content-type"] = "text/csv"
    return output

def run(path_to_pickled_viralcut_collection):
    app.run(host='0.0.0.0', port='8080', debug=True)

if __name__ == '__main__':
    run()