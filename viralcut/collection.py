import os
import json
import pickle
import pandas as pd
from collections import defaultdict
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa

from . import config
from .guide import Guide

class ViralCutCollection:
    def __init__(self, _loading_from_pickled=False):
        self.guides: dict[str, Guide] = {}
        self.target = {}
        self.accessions = []
        self.accession_to_tax_id = {}
        self.node_scores_calculated_for_guides = []
        self.node_scores = {} # key: tuple(node, guide, score_name); value: score
        self.tree_map = {}
        self._ncbi = NCBITaxa()
        self._ncbi_tree = None
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
        self.tree_map = {}
        for node in self._ncbi_tree.traverse("postorder"):
            self.tree_map[str(node.taxid)] = node

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

        for node in tree.traverse("postorder"):
            if node.scored == True:
                continue
            elif len(node.children) > 0:
                for guide in guides:
                    scores = [x.score[guide] for x in node.children]
                    mit_score = [x['mit'] for x in scores if x != -1]
                    cfd_score = [x['cfd'] for x in scores if x != -1]
                    node.score[guide]['mit'] = sum(mit_score) / len(mit_score)
                    node.score[guide]['cfd'] = sum(cfd_score) / len(cfd_score)
                    node.score[guide]['unique_sites'] = sum([x['unique_sites'] for x in scores if x != -1]) 
                    node.score[guide]['total_sites'] = sum([x['total_sites'] for x in scores if x != -1])
            else:
                continue

    def to_pickle(self, filename):
        # Some parts of the Collection cannot be Pickled.
        # They will be removed then regenerated when Unpickled.

        #self._ncbi_tree_newick = write_newick(self._ncbi_tree)

        #del self._ncbi
        #del self._ncbi_tree

        self._pickle_filepath = filename

        with open(filename, 'wb') as fp:
            pickle.dump(self, fp)

    def output_results(self, output_dir):

        target23 = []
        target30 = []
        header = []
        strand = []
        pos = []
        CDE_score = []
        CDE_range = []
        CDE_iqr = []
        CDE_passed = []

        for guide in self.guides.values():
            for i in range(len(guide['30mer'])):
                target23.append(guide.seq)
                target30.append(guide['30mer'][i])
                CDE_score.append(guide['CDE_score'][i])
                CDE_range.append(guide['CDE_UQ_range'][i])
                CDE_iqr.append(guide['CDE_UQ_IQR'][i])
                CDE_passed.append(guide['CDE_passed'][i])
                header.append(guide['header'][i])
                strand.append(guide['strand'][i])
                pos.append(f"{guide['start'][i]}:{guide['end'][i]}")

        ontarget_results = pd.DataFrame({
            'guide': target23,
            '30mer': target30,
            'header': header, 
            'pos': pos,
            'strand': strand,
            'CDE_score': CFDScore, 
            'CDE_range': CFDScore, 
            'CDE_iqr': CFDScore, 
            'CDE_passed': CFDScore
            })

        target23 = []
        accession = []
        CFDScore = []
        MITScore = []
        totalSites = []
        uniqueSites = []

        # Convert tree to dataframe for comparisons
        for guide in self.guides.values():
            for acc, score in guide.assembly_scores.items():
                target23.append(guide.seq)
                accession.append(acc)
                MITScore.append(score['mit'])
                CFDScore.append(score['cfd'])
                totalSites.append(score['total_sites'])
                uniqueSites.append(score['unique_sites'])
        
        offtarget_results = pd.DataFrame({
            'guide': target23, 
            'accession': accession, 
            'MIT score': MITScore, 
            'CFD score': CFDScore, 
            'total_sites': totalSites, 
            'unique_sites': uniqueSites
            })
        offtarget_results['tax_id'] = offtarget_results['accession'].apply(lambda x : self.accession_to_tax_id[x])