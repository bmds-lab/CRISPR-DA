import json
import os
import time
import pickle
from collections import defaultdict
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
from ete3.parser.newick import read_newick, write_newick
import pandas as pd

from . import config
from .guide import Guide

CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNKNOWN = '?'
CODE_ERROR = '!'

class ViralCutCollection:
    def __init__(self, _loading_from_pickled=False):
        self.guides: dict[str, Guide] = {}
        self.target = None
        self.accessions = []
        self.accession_to_tax_id = {}
        self.target_properties = {}
        self.node_scores_calculated_for_guides = []
        self.node_scores = {} # key: tuple(node, guide, score_name); value: score
        self.tree_map = {}
        self._ncbi = NCBITaxa()
        self._ncbi_tree = None
        self._ncbi_tree_newick = None
        self._pickle_filepath = None
        print('ViralCutCollection.__init__ finished')

    def __getitem__(self, key):
        k = self._get_guide_key(key)
        return self.guides[k]

    def __setitem__(self, key, value):
        if key not in self.guides:
            key = self._get_guide_key(key, ignore_nonexist=True)
        self.guides[key] = value

    def __iter__(self):
        for g in self.guides:
            yield g

    def __str__(self):
        props = {
            'num-guides' : len(self.guides)
        }

        if 'geneId' in self.gene_properties:
            props['gene-id'] = self.gene_properties['geneId']

        strProps = ' '.join([
            f"{x}='{props[x]}'"
            for x in props
        ])

        return f'<ViralCutCollection {strProps}>'

    def __getstate__(self):
        # From the Pickle docs:
        # https://docs.python.org/3/library/pickle.html#handling-stateful-objects

        # Copy the object's state from self.__dict__ which contains
        # all our instance attributes. Always use the dict.copy()
        # method to avoid modifying the original state.
        state = self.__dict__.copy()
        
        # Remove the unpicklable entries.
        del state['_ncbi']
        del state['_ncbi_tree']
        return state

    def __setstate__(self, state):
        # Restore instance attributes (i.e., filename and lineno).
        self.__dict__.update(state)

        self._ncbi = NCBITaxa()
        self._ncbi_tree = None
        self._ncbi_tree = self.get_ncbi_tax_tree()
        
        

    def _autoSaveIfPossible(self):
        did_auto_save = False

        if self._pickle_filepath is not None:
            if os.path.exists(self._pickle_filepath):
                if config.VERBOSE:
                    print(f'Autosaving to: {self._pickle_filepath}')

                try:
                    self.to_pickle(self._pickle_filepath)
                    did_auto_save = True
                except Exception as e:
                    raise e

        if not did_auto_save and config.VERBOSE:
            print('Could not autosave ViralCutCollection')

    def _get_guide_key(self, key, ignore_nonexist=False):
        if key in self.guides:
            return key

        # an exact match did not exist. check for substrings.
        # usually this is helpful if a PAM-less sequence is provided
        # i.e., ATCGATCGATCGATCGATCGAGG versus ATCGATCGATCGATCGATCG
        keys = [x for x in self.guides if x[0:len(key)] == key]

        if len(keys) > 1:
            if not ignore_nonexist:
                raise RuntimeError(f'Multiple matches to {key} in collection')
        elif len(keys) == 1:
            return keys[0]
        else:
            return key

    def update_phylogenetic_tree(self, tax_ids):
        ''' 
        Updates the phylogenetic tree to the smallest tree that connects all your query taxids
        Returns:
            None
        '''
        self._ncbi_tree = self._ncbi.get_topology(tax_ids, intermediate_nodes=False)
        self.map_phylogentic_tree()


    def reset_phylogenetic_tree(self, root_tax_id):
        '''
        Resets the phylogenetic tree using the root tax id
        '''
        self._ncbi_tree = self._ncbi.get_topology(root_tax_id, intermediate_nodes=True)
        self.map_phylogentic_tree()

    def map_phylogentic_tree(self):
        self.tree_map = {}
        for node in self._ncbi_tree.traverse("postorder"):
            self.tree_map[str(node.taxid)] = node


    def get_tax_ids_in_analysis(self):
        '''Fetches a list of NCBI taxonomy IDs that exist in the Collection

        Returns:
            A list of integers
        '''

        return list(map(int, self.accession_to_tax_id.values()))


    def guides_to_dataframe(self):
        data = {'seq' : []}
        all_props = set()

        for g in self.guides:
            for prop in self.guides[g].props:
                all_props.add(prop)
            
            data['seq'].append(g)
        
            for prop in all_props:
                if prop not in data:
                    data[prop] = []
        
                if prop in self.guides[g].props:
                    data[prop].append(self.guides[g].props[prop])
                else:
                    data[prop].append(CODE_UNKNOWN)
        
        return pd.DataFrame(data)

    def assembly_scores_to_dataframe(self):
        '''This function takes the assembly scores of each guide (dict of lists) and merges them
        into one data structure (dict of lists). Adding the extract column for the guide identity is
        needed. The DataFrame will look something like:

        > guide        accession score_name       score unique_sites total_sites
        > 0  TACACTAATTCTTTCACACGTGG  GCA_000859285.1        mit  100.000000     0.000000    0.000000
        > 1  TACACTAATTCTTTCACACGTGG  GCA_000886535.1        mit  100.000000     0.000000    0.000000
        > 2  TACACTAATTCTTTCACACGTGG  GCA_000884175.1        cfd  100.000000     0.000000    0.000000


        Returns:
            A DataFrame
        '''

        scores = {'guide' : []}

        # for guide in self.guides:

        #     for k in self.guides[guide].assembly_scores:
        #         if k not in scores:
        #             scores[k] = []

        #         scores[k] += self.guides[guide].assembly_scores[k]

        #     scores['guide'] += [guide] * len(self.guides[guide].assembly_scores[k])

        return pd.DataFrame(scores)

    def calculate_node_scores(self):
        '''
        After completing the off-target scoring, this method will add the scores to the tree.
        It will then propergate the scores (by averaging it's children) through the tree.
        '''
        tree = self._ncbi_tree

        if config.VERBOSE:
            print(f'Preparing scores data structure')

        guides = [y for y, x in self.guides.items() if x.assembly_scores != {}]

        if config.VERBOSE:
            n = len(list(tree.traverse(strategy="postorder")))
            print(f'Calculating scores for {n} nodes using depth-first traversal')

        for node in tree.traverse("postorder"):
            node.score = defaultdict(lambda : {'mit': -1, 'cfd': -1, 'unique_sites': -1, 'total_sites': -1})
            node.scored = False

        for accession, taxid in self.accession_to_tax_id.items():
            node = self.tree_map[str(taxid)]
            for guide in guides:
                try:
                    mit, cfd, unique_sites, total_sites = self.guides[guide].assembly_scores[accession].values()
                    node.score[guide]['mit'] = float(mit)
                    node.score[guide]['cfd'] = float(cfd)
                    node.score[guide]['unique_sites'] = int(unique_sites)
                    node.score[guide]['total_sites'] = int(total_sites)
                    node.scored = True
                except:
                    continue

        guide = guides[0]
        for node in tree.traverse("postorder"):
            if node.scored == True:
                continue
            elif len(node.children) > 0:
                scores = [x.score[guide] for x in node.children]
                mit_score = [x['mit'] for x in scores if x != -1]
                cfd_score = [x['mit'] for x in scores if x != -1]
                node.score[guide]['mit'] = sum(mit_score) / len(mit_score)
                node.score[guide]['cfd'] = sum(cfd_score) / len(cfd_score)
                node.score[guide]['unique_sites'] = sum([x['unique_sites'] for x in scores if x != -1]) 
                node.score[guide]['total_sites'] = sum([x['total_sites'] for x in scores if x != -1])
            else:
                continue

    def get_node_score(self, tax_id, guide, score_name, r=False):
        # have the scores for this guide been calculated?
        # if not, calculate them now

        key = (tax_id, guide, score_name)
        #print(key)
        try:
            s = 10_000 / (100.0 + self.node_scores[key])
            #print('got score')
            print(key, s)
            return s
        except KeyError as e:
            print('key error when getting node score')
            pass

        scores_for_guide_have_been_calculated = (score_name, guide) in self.node_scores_calculated_for_guides

        if r:
            # already been recursively called. if we have reached here,
            # then the score really doesn't exist.
            print('get_node_score really couldn\'t find a score')
            return 0

        if not scores_for_guide_have_been_calculated:
            print('get_node_score is calling for scores to be calculated')
            self.calculate_node_scores(guides=[guide])
            return self.get_node_score(tax_id, guide, score_name, r=True)
        
        return -1


    def get_descendant_tax_ids_from_root_tax_id(self, tax_id=None, max_depth=None, max_nodes=None, include_level_zero=True):
        '''Given some taxonomy ID, find the taxonomy IDs of descendant nodes

        Arguments:
            root_tax_id (int): The root tax_id. Optional. default is config.ROOT_TAX_ID

        Returns:
            A list of tax IDS

        '''

        tree = self._ncbi_tree
        include_nodes = set()
        depth_offset = 0

        found_tax_id_in_tree = False
        if tax_id is not None:
            include_nodes.update([tax_id])

            # a tax_id has been specified, first find it
            for idx, i in enumerate(tree.traverse(strategy="levelorder")):
                if int(i.name) == tax_id:
                    found_tax_id_in_tree = True
                    print(f'updating root of tree from {tree.name} to {i.name}')
                    tree = i
                    depth_offset = len(i.get_ancestors()) - 1
                    break

            if not      found_tax_id_in_tree:
                print(f'Could not find tax_id {tax_id} in the tree')
                return [tax_id]

        for idx, i in enumerate(tree.traverse(strategy="levelorder")):
            depth = len(i.get_ancestors()) - 1 - depth_offset

            doBreak = False
            doBreak |= (max_nodes is not None and len(include_nodes) > max_nodes)
            doBreak |= (max_depth is not None and depth > max_depth)

            if doBreak:
                break

            doAdd = i.up is not None
            doAdd &= (
                (include_level_zero and depth >= 0) or
                (not include_level_zero and depth > 0)
            )
            if doAdd:
                #print(f'depth is {depth}')
                include_nodes.update([
                    int(node.name)
                    for node in i.up.children
                ])

        return list(include_nodes)

    def generate_newick_string_from_tax_ids(self, tax_ids):
        '''This function generates Newick string of the data needed to visualise the
        phylogenetic tree in a web browser. See https://en.wikipedia.org/wiki/Newick_format

        Arguments:
            tax_ids (list):  The taxonomy IDs to include in the tree

        Returns:
            A Newick string
        '''

        return write_newick(
            self._ncbi.get_topology(
                tax_ids
            )
        )

    def to_pickle(self, filename):
        # Some parts of the Collection cannot be Pickled.
        # They will be removed then regenerated when Unpickled.

        #self._ncbi_tree_newick = write_newick(self._ncbi_tree)

        #del self._ncbi
        #del self._ncbi_tree

        self._pickle_filepath = filename

        with open(filename, 'wb') as fp:
            pickle.dump(self, fp)

    def to_human_readable_files(self, filename):
        df_guides = self.guides_to_dataframe()
        
        df_assembly_scores = self.assembly_scores_to_dataframe()
        
        df_node_scores = pd.DataFrame(self.node_scores).T

        with open(f"{filename}-guides.csv", "w") as fp:
            df_guides.to_csv(fp, index=False)

        with open(f"{filename}-guides.md", "w") as fp:
            df_guides.to_markdown(fp, index=False)

        with open(f"{filename}-gene-properties.txt", "w") as fp:
            fp.write(
                json.dumps(
                    self.gene_properties,
                    sort_keys=True,
                    indent=4
                )
            )

        with open(f"{filename}-assembly-scores.csv", "w") as fp:
            df_assembly_scores.to_csv(fp, index=False)

        with open(f"{filename}-node-scores.csv", "w") as fp:
            df_node_scores.to_csv(fp, index=False)