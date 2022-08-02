from . import analysis
from . import data
from .data import collection_from_pickle

import pandas as pd
from flask import Flask, send_from_directory, make_response, jsonify

app = Flask(__name__,
    static_url_path='', 
    static_folder='viewer',
    template_folder='templates'
)

MAX_TAX_IDS = 300
VC_COLLECTION = None

#accessions = data.get_accessions_from_ncbi_table_export('/home/jake/ViralCut/examples/viruses.csv')[0:50]
#tax_ids = analysis.get_tax_ids_from_accessions(accessions) + 

def get_tax_ids(root_tax_id):
    #return data.get_cached_tax_ids()[:MAX_TAX_IDS]
    tax_ids = VC_COLLECTION.get_descendant_tax_ids_from_root_tax_id(
        tax_id=root_tax_id, 
        max_depth=3, 
        max_nodes=MAX_TAX_IDS
    )
    return tax_ids[:MAX_TAX_IDS]
    

@app.route("/")
def root():
    return send_from_directory('viewer', 'index.html')

@app.route("/newick/<int:root_tax_id>")
def newick(root_tax_id):
    tax_ids = get_tax_ids(root_tax_id)
    print('tax_ids: ', tax_ids)
    output = make_response(
        analysis.generate_newick_string_from_tax_ids(tax_ids)
    )
    output.headers["Content-type"] = "text/plain"
    return output
    
@app.route("/subnewick/<int:root_tax_id>")
def subnewick(root_tax_id):
    tax_ids = VC_COLLECTION.get_descendant_tax_ids_from_root_tax_id(
        tax_id=root_tax_id, 
        max_depth=4, 
        max_nodes=MAX_TAX_IDS,
        include_level_zero=False,
    )
    
    output = make_response(
        analysis.generate_newick_string_from_tax_ids(tax_ids)
    )
    output.headers["Content-type"] = "text/plain"
    return output
    

@app.route("/scores/<int:root_tax_id>/<guide>")
def scores(root_tax_id, guide):
    #tax_ids = get_tax_ids(root_tax_id)
    #df = analysis.generate_df_phylo_node_scores_from_tax_ids(tax_ids)
    
    #df = VC_COLLECTION.calculate_node_scores(
    #    root_tax_id=root_tax_id,
    #    guides=[guide]
    #)
    
    #tax_ids = get_tax_ids(root_tax_id)
    #                                                                 
    #for tax_id in tax_ids:
    #    print(tax_id)
    #    VC_COLLECTION.get_node_score(tax_id, guide, 'mit')
    #
    #df = pd.concat([
    #    VC_COLLECTION.get_node_score(tax_id, guide, 'mit')
    #    for tax_id in tax_ids
    #])

    df = VC_COLLECTION.node_scores

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
    
@app.route("/subspecies/<int:root_tax_id>")
def subspecies(root_tax_id):
    tax_ids = VC_COLLECTION.get_descendant_tax_ids_from_root_tax_id(
        tax_id=root_tax_id, 
        max_depth=4, 
        max_nodes=MAX_TAX_IDS,
        include_level_zero=False,
    )
    print('fetching subspecies annotations: ', tax_ids)
    df = analysis.generate_df_of_species_data_from_tax_ids(tax_ids)
    output = make_response(
        df.to_csv(index=False)
    )
    
    output.headers["Content-type"] = "text/csv"
    return output
    
@app.route("/guides")
def guides():
    df = pd.DataFrame(VC_COLLECTION.guides.keys(), columns=['guide'])
    output = make_response(
        df.to_csv(index=False)
    )
    output.headers["Content-type"] = "text/csv"
    return output

@app.route("/subtree/<int:root_tax_id>/<guide>/<score>")
def subtree(root_tax_id, guide, score):
    returns = {
        'root_tax_id' : root_tax_id,
        'newick' : None,
        'annotations' : None,
        'scores' : {},
        'leaves_with_children' : {}
    }
    
    tax_ids = VC_COLLECTION.get_descendant_tax_ids_from_root_tax_id(
        tax_id=root_tax_id, 
        max_depth=2, 
        max_nodes=MAX_TAX_IDS,
        include_level_zero=False,
    )
    
    # generate newick string
    returns['newick'] = analysis.generate_newick_string_from_tax_ids(tax_ids)
    
    # generate annotation info
    returns['annotations'] = analysis.generate_df_of_species_data_from_tax_ids(tax_ids).to_dict()
    
    # which leaves can be expanded further?
    
    
    # scores
    import random
    
    VC_COLLECTION.calculate_node_scores(
        guides=[guide],
        root_tax_id=root_tax_id
    )
    
    print('getting node scores')
    for tax_id in tax_ids:
        try:
            returns['scores'][tax_id] = VC_COLLECTION.get_node_score(tax_id, guide, score)
            #returns['scores'][tax_id] = random.random() * 100 #VC_COLLECTION.get_node_score(tax_id, guide, 'mit')
        except Exception as e:
            raise e
            returns['scores'][tax_id] = 0
            
    return returns

def run(
    path_to_pickled_viralcut_collection, 
    host='0.0.0.0', 
    port=8080
):
    global VC_COLLECTION
    VC_COLLECTION = collection_from_pickle(path_to_pickled_viralcut_collection)
    VC_COLLECTION.node_scores_calculated_for_guides = []
    VC_COLLECTION.node_scores = {}
    app.run(host=host, port=port, debug=True)

if __name__ == '__main__':
    print('This is not a script.\nYou need to run the server via the ViralCut API or\n   the CLI command that was installed with ViralCut.')
    