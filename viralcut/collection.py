import json
import os
import time
import pickle
from collections import defaultdict
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
from ete3.parser.newick import read_newick, write_newick
import pandas as pd

from . import config

CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNKNOWN = '?'
CODE_ERROR = '!'

class ViralCutCollection:
    def __init__(self, _loading_from_pickled=False):
        self.guides = {}
        self.target = None
        self.accessions = []
        self.accession_to_tax_id = {}
        self.gene_properties = {}
        self.node_scores_calculated_for_guides = []
        self.node_scores = {} # key: tuple(node, guide, score_name); value: score
        self._ncbi = None
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

    def reset_phylogenetic_tree(self, root_tax_id):
        '''
        Resets the phylogenetic tree using the root tax id
        '''
        self._ncbi = NCBITaxa()
        self._ncbi_tree = self._ncbi.get_topology(root_tax_id, intermediate_nodes=True)


    def get_tax_ids_in_analysis(self):
        '''Fetches a list of NCBI taxonomy IDs that exist in the Collection

        Returns:
            A list of integers
        '''

        return list(map(int, self.accession_to_tax_id.values()))

    def get_ncbi_tax_tree(self):
        if self._ncbi_tree is None:
            self._ncbi_tree = self._ncbi.get_topology(
                self.get_tax_ids_in_analysis(),
                intermediate_nodes=True
            )
        return self._ncbi_tree

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

    def calculate_node_scores(self,
        guides=None,
        root_tax_id=None,
        autosave=True
    ):
        '''
        Arguments:
            autosave (bool): autosave but only if new scores are calculated.
        
        
        '''
        if root_tax_id is None:
            root_tax_id = config.ROOT_TAX_ID

        if config.VERBOSE:
            print(f'Calculating node scores, with root taxonomy ID: {root_tax_id}')

        # tree = self.get_ncbi_tax_tree()
        tree = self._ncbi_tree
        if config.VERBOSE:
            print(f'Preparing scores data structure')


        # df_accs_to_tax_id = pd.DataFrame(self.accession_to_tax_id.items(), columns=['accession', 'tax_id'])
        # '''
        # df_accs_to_tax_id is structured as following:
        #    |    | accession       |   tax_id |
        #    |---:|:----------------|---------:|
        #    |  0 | GCA_000859285.1 |    85708 |
        #    |  1 | GCA_000886535.1 |    11676 |
        #    |  2 | GCA_000884175.1 |  1645793 |
        # '''


        # df_scores = self.assembly_scores_to_dataframe()
        # '''
        # df_scores is structured as following:
        #    |    | guide                   | accession       | score_name   | score | unique_sites | total_sites |  tax_id |
        #    |---:|:------------------------|:----------------|:-------------|------:|-------------:|------------:|--------:|
        #    |  0 | TACACTAATTCTTTCACACGTGG | GCA_000859285.1 | mit          |   100 |            0 |           0 |   85708 |
        #    |  1 | TACACTAATTCTTTCACACGTGG | GCA_000886535.1 | mit          |   100 |            0 |           0 |   11676 |
        #    |  2 | TACACTAATTCTTTCACACGTGG | GCA_000884175.1 | cfd          |   100 |            0 |           0 | 1645793 |
        # '''

        # df_accs_to_tax_id.set_index('accession')

        # possible_guides = list(set(df_scores['guide']))

        # df_scores.set_index(['accession', 'guide'])

        # df_scores = df_scores.merge(df_accs_to_tax_id, on='accession', how='left')
        # score_names = ['mit'] #set(df_scores['score_name']):

        # df_scores.set_index(['tax_id', 'score_name', 'guide'], inplace=True)
        # df_scores.sort_index(inplace=True)

        # # df_scores now looks like:
        # '''
        # | (tax_id, score_name, guide)               | accession       |   score |   unique_sites |   total_sites |
        # |:------------------------------------------|:----------------|--------:|---------------:|--------------:|
        # | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_018595055.1 |       0 |              0 |             0 |
        # | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_900323395.1 |       0 |              0 |             0 |
        # | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_006458985.1 |       0 |              0 |             0 |
        # '''

        # print('Printing rows with non-zero scores')
        # print(df_scores[df_scores['score'] > 0])

        # if guides is None:
        #     guides = possible_guides#[possible_guides[1]]

        guides = [y for y, x in self.guides.items() if x.assembly_scores != {}]

        if config.VERBOSE:
            n = len(list(tree.traverse(strategy="postorder")))
            print(f'Calculating scores for {n} nodes using depth-first traversal')

        mapped_tree = {}

        # depth_to_node = defaultdict(list)
        for node in tree.traverse("postorder"):
            # node.score = defaultdict(lambda : {'mit': -1, 'cfd': -1})
            node.score = defaultdict(lambda : {'mit': [], 'cfd': []})
            node.scored = False
            mapped_tree[node.taxid] = node
        #     dist = tree.get_distance(node)
        #     depth_to_node[dist].append(node)
        # depths = [x for x in depth_to_node.keys()]
        # depths.sort(reverse=True)

        start_time = time.time()

        for accession, taxid in self.accession_to_tax_id.items():
            node = mapped_tree[taxid]
            for guide in guides:
                try:
                    guide_results = self.guides[guide].assembly_scores[accession]
                    score = guide_results['score']
                    # node.score[guide][guide_results['score_name']] = score
                    node.score[guide][guide_results['score_name']].append((accession, score))
                    # if node.scored:
                    #     print(f'{taxid} has been scored multiple times')
                    node.scored = True
                    # propagate change up the tree to the root
                    while True:
                        node = node.up
                        if node is None:
                            break
                        # node.score[guide][guide_results['score_name']] = max(node.score[guide][guide_results['score_name']], score)
                        node.score[guide][guide_results['score_name']].append((accession, score))
                except:
                    continue



        # # 
        # for d in depths:
        #     print(d)
        #     for node in depth_to_node[d]:
        #         # Leaf
        #         if len(node.children) == 0:
        #             continue
        #         # propergate scores up from children
        #         else:
        #             for guide in guides:
        #                 node.score[guide]['mit'] = max([x.score[guide]['mit'] for x in node.children])



        # start_time = time.time() 3925.530910253525

        # for d in depths:
        #     print(d)
        #     for node in depth_to_node[d]:
        #         node.score = defaultdict(lambda : {'mit': [], 'cfd': []})
        #         # Leaf, get score from dataframe
        #         if len(node.children) == 0:
        #             for guide in guides:
        #                 for score_name in config.ISSL_SCORES:
        #                     try:
        #                         issl_score = df_scores.loc[(node.taxid, score_name, guide), 'score']
        #                         node.score[guide][score_name].append(issl_score)
        #                     except KeyError as e:
        #                         node.score[guide][score_name].append('?')
        #         # Get score from children
        #         else:
        #             for guide in guides:
        #                 for score_name in config.ISSL_SCORES:
        #                     for child in node.children:
        #                         for s in child.score[guide][score_name]:
        #                             node.score[guide][score_name].append(s)

        # start_time = time.time()

        # # depth-first traversal 10025.991832017899
        # new_scores_calculated = False
        # for score_name in score_names:

        #     for guide in guides:
        #         key = (score_name, guide)
        #         if key in self.node_scores_calculated_for_guides:
        #             print(f'Scores already calculated for: {key}')
        #             continue
        #         else:
        #             self.node_scores_calculated_for_guides.append(
        #                 key
        #             )
        #             new_scores_calculated = True

        #         for idx, i in enumerate(tree.traverse(strategy="postorder")):

        #             if i.name == '':
        #                 continue

        #             tax_id = int(i.name)

        #             if len(i.children) == 0:
        #                 # is a leaf so find minimum of assembly scores

        #                 try:
        #                     scores = df_scores.loc[(tax_id, score_name, guide), 'score']
        #                     score = scores.sum()
        #                 except KeyError as e:
        #                     score = -101
                            
                            
        #                 if score > 0:
        #                         print(tax_id, score_name, guide, '\n', scores[scores > 0], '\n\n\n')
                        
        #                 # propagate change up the tree to the root
        #                 node = i.up
        #                 while True:
        #                     if node is None:
        #                         break
                                
        #                     if hasattr(node, 'score'):
        #                         node.score = max(
        #                             node.score,
        #                             score
        #                         )
        #                     else:
        #                         node.score = score

        #                     #print((
        #                     #    int(node.name),
        #                     #    guide,
        #                     #    score_name
        #                     #), score)
        #                     self.node_scores[(
        #                         int(node.name),
        #                         guide,
        #                         score_name
        #                     )] = score
                            
        #                     node = node.up
                            
        #                     if node is None:
        #                         break

        # if autosave and new_scores_calculated:
        #     self._autoSaveIfPossible()
        
        elapsed_time = time.time() - start_time
        print(elapsed_time)
        return True

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