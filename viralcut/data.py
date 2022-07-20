import csv
import json
import os
import subprocess
import multiprocessing
import math
import zipfile as zf
from glob import glob
from collections import deque
from io import BytesIO, TextIOWrapper

try:
   import cPickle as pickle
except:
   import pickle
   
from . import config

import ncbi.datasets
from ete3.ncbi_taxonomy.ncbiquery import NCBITaxa
from ete3.parser.newick import read_newick, write_newick
import pandas as pd
import tqdm

class Guide:
    def __init__(self, seq):
        self.seq = seq

        # properties of a guide
        self.props = {}
        
        # scores
        self.assembly_scores = {
            'accession' : [],
            'score_name' : [],
            'score' : [],
            'unique_sites' : [],
            'total_sites' : []
        }

    def __getitem__(self, key):
        if key not in self.props:
            self.props[key] = CODE_UNKNOWN
        return self.props[key]

    def __setitem__(self, key, value):
        self.props[key] = value

    def __str__(self):
        strProps = ' '.join([
            f"{x}='{self.props[x]}'"
            for x in self.props
        ])
        return f"<Guide seq='{self.seq}' {strProps}>"
        
    def add_assembly_score(self, accession, score_name, score, unique_sites, total_sites):
        self.assembly_scores['accession'].append(accession)
        self.assembly_scores['score_name'].append(score_name)
        self.assembly_scores['score'].append(score)
        self.assembly_scores['unique_sites'].append(unique_sites)
        self.assembly_scores['total_sites'].append(total_sites)
        
    def assembly_scores_to_dataframe(self):
        return pd.DataFrame(self.assembly_scores.items())

class ViralCutCollection:
    def __init__(self):
        self.guides = {}
        self.gene_id = None
        self.accessions = []
        self.accession_to_tax_id = {}
        self.gene_properties = {}
        self.node_scores = pd.DataFrame()
        self._ncbi = None
        self._ncbi_tree = None
        self._ncbi_tree_newick = None

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
        dfs = []        
        for g in self.guides:
            df = g.assembly_scores_to_dataframe()
            df['guide'] = g
            dfs.append(df)
        return pd.concat(dfs)
    
        #data = {'seq' : []}
        #       all_props = set()
        #       
        #       for g in self.guides:
        #           for prop in self.guides[g].props:
        #               all_props.add(prop)
        #     data['seq'].append(g)
        #
        #    for prop in all_props:
        #        if prop not in data:
        #            data[prop] = []
        #            
        #        if prop in self.guides[g].props:
        #            data[prop].append(self.guides[g].props[prop])
        #        else:
        #            data[prop].append(CODE_UNKNOWN)
        #
        #return pd.DataFrame(data)

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
        
        for guide in self.guides:
        
            for k in self.guides[guide].assembly_scores:
                if k not in scores:
                    scores[k] = []
                    
                scores[k] += self.guides[guide].assembly_scores[k]
            
            scores['guide'] += [guide] * len(self.guides[guide].assembly_scores[k])
            
        return pd.DataFrame(scores)

    def calculate_node_scores(self,
        guides=None,
        root_tax_id=None
    ):
        if root_tax_id is None:
            root_tax_id = config.ROOT_TAX_ID

        if config.VERBOSE:
            print(f'Preparing tree, with root: {root_tax_id}')

        tree = self.get_ncbi_tax_tree()
        
        if config.VERBOSE:
            print(f'Preparing scores data structure')
        
        df_scores = self.assembly_scores_to_dataframe()
        
        '''
        The following DataFrame is structured as following:
           |    | accession       |   tax_id |
           |---:|:----------------|---------:|
           |  0 | GCA_000859285.1 |    85708 |
           |  1 | GCA_000886535.1 |    11676 |
           |  2 | GCA_000884175.1 |  1645793 |
        '''
        df_accs_to_tax_id = pd.DataFrame(self.accession_to_tax_id.items(), columns=['accession', 'tax_id'])
        
        '''
        The following DataFrame is structured as following:
           |    | guide                   | accession       | score_name   | score | unique_sites | total_sites |  tax_id |
           |---:|:------------------------|:----------------|:-------------|------:|-------------:|------------:|--------:|
           |  0 | TACACTAATTCTTTCACACGTGG | GCA_000859285.1 | mit          |   100 |            0 |           0 |   85708 |
           |  1 | TACACTAATTCTTTCACACGTGG | GCA_000886535.1 | mit          |   100 |            0 |           0 |   11676 |
           |  2 | TACACTAATTCTTTCACACGTGG | GCA_000884175.1 | cfd          |   100 |            0 |           0 | 1645793 |
        '''
        
        df_accs_to_tax_id.set_index('accession')
        df_scores.set_index(['accession', 'guide'])
        df_scores = df_scores.merge(df_accs_to_tax_id, on='accession', how='left')
        
        score_names = ['mit'] #set(df_scores['score_name']):

        df_scores.set_index(['tax_id', 'score_name', 'guide'], inplace=True)
        df_scores.sort_index(inplace=True)
    
        # df_scores now looks like:
        '''
        | (tax_id, score_name, guide)               | accession       |   score |   unique_sites |   total_sites |
        |:------------------------------------------|:----------------|--------:|---------------:|--------------:|
        | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_018595055.1 |       0 |              0 |             0 |
        | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_900323395.1 |       0 |              0 |             0 |
        | (10243, 'cfd', 'AAACCTAGCCAAATGTACCATGG') | GCA_006458985.1 |       0 |              0 |             0 |
        '''
        
        df_scores.to_csv('temp-20220719.txt', sep='\t')
        if guides is None:
            guides = set(df_scores['guide'])
            guides = [list(guides)[0]]
        
        if config.VERBOSE:
            print(f'Indexing scores')

        data = {'tax_id': [], 'score_name' : [], 'guide' : [], 'score' : []}
        
        if config.VERBOSE:
            print(f'Calculating node scores using depth-first traversal')
        
        # depth-first traversal
        for score_name in score_names:
            
            for guide in guides:
                
                tax_id_to_score = {}
                for idx, i in enumerate(tree.traverse(strategy="postorder")):
                    tax_id = int(i.name)
                    
                    if len(i.children) == 0:
                        # is a leaf so find minimum of assembly scores

                        try:
                            scores = df_scores.loc[(tax_id, score_name, guide), 'score']
                            #print(tax_id, score_name, guide, '\n', scores, '\n\n\n')
                            score = scores.max()
                        except KeyError as e:
                            score = -101
                            
                                    
                    else:
                        
                        # find the minimum of child node scores
                        immediate_child_tax_ids = [int(child.name) for child in i.children]

                        if len(immediate_child_tax_ids) == 0:
                            score = math.nan
                        else:
                            score = max([
                                tax_id_to_score[tax_id]
                                for tax_id in immediate_child_tax_ids
                                if not math.isnan(tax_id_to_score[tax_id])
                            ], default=math.nan)
                            
                        score = 10_000.0 / (100.0 + score)

                    tax_id_to_score[tax_id] = score
                    data['score_name'].append(score_name)
                    data['tax_id'].append(tax_id)
                    data['guide'].append(guide)
                    data['score'].append(score)
            
                    try:
                        print('\t', score_name, tax_id, i.up.name, guide, score)
                    except Exception as e:
                        pass
                        
        if len(self.node_scores) == 0:
            self.node_scores = pd.DataFrame(data)
        else:
            self.node_scores = pd.concat([
                self.node_scores, 
                pd.DataFrame(data)
            ])
        
        print(self.node_scores)
    
    def get_node_score(self, tax_id, guide, score_name):
        if (tax_id, guide, score_name) not in self.node_scores:
            self.calculate_node_scores(guides=[guide])
        return self.node_scores[(tax_id, guide, score_name)]

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
        
        self._ncbi_tree_newick = write_newick(self._ncbi_tree)
        
        del self._ncbi
        del self._ncbi_tree
    
        with open(filename, 'wb') as fp:
            pickle.dump(self.__dict__, fp)
    
    def to_human_readable_files(self, filename):
        df_guides = self.guides_to_dataframe()
        df_assembly_scores = self.assembly_scores_to_dataframe()
    
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
 
def collection_from_pickle(filename):
    collection = ViralCutCollection()
    with open(filename, 'rb') as fp:
        collection.__dict__ = pickle.load(fp)
    
    # Recreate parts of the Collection that would not have been Pickled
    collection._ncbi = NCBITaxa()
    
    try:
        collection._ncbi_tree = read_newick(collection._ncbi_tree_newick)
        raise RuntimeError('test')
    except Exception as e:
        collection._ncbi_tree = collection._ncbi.get_topology(
            collection.get_tax_ids_in_analysis(),
            intermediate_nodes=True
        )
    
    return collection

def parse_fna(stream):
    '''Parse some iterable object as a multi-FASTA file.
    Yield after reading each FASTA block.

    Arguments:
        stream (iterable):  An iterable object to read

    '''
    header = None
    seqs = []
    for line in stream:
        line = line.strip()

        if line[0] == '>':
            if header is not None:
                yield header, ''.join(seqs)
            header = line
        else:
            seqs.append(line)
    yield header, ''.join(seqs)


def download_ncbi_genes(gene_ids):
    '''Downloads the sequence of the specified gene

    Arguments:
        gene_ids (list): List of NCBI gene IDs as integers

    Returns:
        A list of downloaded genes by ID
    '''

    # Convert gene_ids to list if non-list type is provided
    if not isinstance(gene_ids, list):
        gene_ids = [int(gene_ids)]

    # Some genes may already be cached,
    genes_to_download = [
        gene_id
        for gene_id in gene_ids
        if not os.path.exists(get_gene_cache(gene_id, mkdir=False))
    ]

    # which means we may not have anything to download.
    if not len(genes_to_download):
        if config.VERBOSE:
            print('No genes to download')
        return []

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(gene_ids) != len(genes_to_download):
        print((
            f'{len(gene_ids) - len(genes_to_download)} of {len(gene_ids)} requested exist in the '
            f'cache already.'
        ))

    # Keep track of which genes were actually downloaded.
    gene_ids_downloaded = []

    # New instance of the NCBI API interface.
    api_instance = ncbi.datasets.GeneApi(ncbi.datasets.ApiClient())

    api_response = api_instance.download_gene_package(
        [int(x) for x in gene_ids],
        include_annotation_type=['FASTA_GENE'],
        _preload_content=False
    )

    # Process the downloaded data without writing it to disk, yet.
    in_memory = BytesIO(api_response.data)
    with zf.ZipFile(in_memory) as zfp:
        files = zfp.namelist()

        if 'ncbi_dataset/data/gene.fna' not in files:
            raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/gene.fna`')

        # Extract the data report
        if 'ncbi_dataset/data/data_report.jsonl' not in files:
            raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/data_report.jsonl`')

        data_report_by_gene_id = {}
        
        with zfp.open('ncbi_dataset/data/data_report.jsonl') as fp:
            for line in fp:
                report = json.loads(line.strip())
                id = report['geneId']
                data_report_by_gene_id[id] = report

        # ZipFile.open() reads as a binary stream but we need a text stream.
        # https://docs.python.org/3/library/io.html#io.TextIOWrapper
        with TextIOWrapper(zfp.open('ncbi_dataset/data/gene.fna')) as fp:
            for header, seq in parse_fna(fp):
                header_gene_id = get_properties_from_ncbi_fasta_header(header, key='GeneID')

                filename = f"{header_gene_id}.fna"
                cache_dir = get_gene_cache(header_gene_id)
                cached_file = os.path.join(cache_dir, filename)
                cache_data_report = os.path.join(cache_dir, 'data_report.json')

                with open(cached_file, 'w') as fp:
                    fp.write(f'{header}\n')
                    fp.write(f'{seq}\n')

                with open(cache_data_report, 'w') as fpW:
                    json.dump(data_report_by_gene_id[header_gene_id], fpW)
                
                if config.VERBOSE:
                    print(f'Downloaded to: {cache_dir}')

                gene_ids_downloaded.append(header_gene_id)

    return gene_ids_downloaded


def download_ncbi_assemblies(accessions, keep_exts=['fna'], merge=False):
    '''Download files associated with the provided accessions to a cache.

    Arguments:
        accessions (list): List of NCBI accessions as strings
        keep_exts (list):  Files with these extensions will be cached
        merge (bool):      Only cache files that have not been cached already

    Return:
        A list of accessions actually downloaded

    '''

    # Convert gene_ids to list if a string is provided
    if isinstance(accessions, str):
        accessions = [accessions]

    # Some accessions may already be cached,
    accs_to_download = [
        accs
        for accs in accessions
        if not os.path.exists(get_assembly_cache(accs, mkdir=False)) or merge
    ]

    # which means we may not have anything to download.
    if not len(accs_to_download):
        if config.VERBOSE:
            print('No assemblies to download')
        return []

    # Notify if only some are being downloaded.
    if config.VERBOSE and len(accessions) != len(accs_to_download):
        print((
            f'{len(accessions) - len(accs_to_download)} of {len(accessions)} '
            f'requested exist in the cache already.'
        ))

    # Keep track of which accessions were actually downloaded.
    accessions_downloaded = []

    # New instance of the NCBI API interface.
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

    for start in range(0, len(accs_to_download), config.NCBI_ACCESSION_BATCH_SIZE):
        api_response = api_instance.download_assembly_package(
            accs_to_download[start:start+config.NCBI_ACCESSION_BATCH_SIZE],
            _preload_content=False
        )

        # Process the downloaded data without writing it to disk, yet.
        in_memory = BytesIO(api_response.data)
        with zf.ZipFile(in_memory) as zfp:
            files = zfp.namelist()
            
            # Extract the data report
            if 'ncbi_dataset/data/assembly_data_report.jsonl' not in files:
                raise RuntimeError('Critical file not provided by NCBI: `ncbi_dataset/data/assembly_data_report.jsonl`')

            data_report_by_accession = {}
            
            with zfp.open('ncbi_dataset/data/assembly_data_report.jsonl') as fp:
                for line in fp:
                    report = json.loads(line.strip())
                    accs = report['assemblyInfo']['assemblyAccession']
                    data_report_by_accession[accs] = report

            for path in files:

                # Be selective in which files get cached.
                accs = os.path.basename(os.path.dirname(path))
                filename = os.path.basename(path)
                cache_dir = get_assembly_cache(accs)
                cached_file = os.path.join(cache_dir, filename)
                cache_data_report = os.path.join(cache_dir, 'data_report.json')

                if not('GCF' in accs or 'GCA' in accs):
                    continue

                if filename.split('.')[-1] in keep_exts:
                
                    if not os.path.exists(cached_file):
                        with zfp.open(path) as fp, open(cached_file, 'wb') as fpW:
                            fpW.writelines(fp.readlines())
                        
                        if config.VERBOSE:
                            print(f'Downloaded to: {cached_file}')
                    
                    if not os.path.exists(cache_data_report):
                        with open(cache_data_report, 'w') as fpW:
                            json.dump(data_report_by_accession[accs], fpW)
                        
                        if config.VERBOSE:
                            print(f'Downloaded to: {cache_data_report}')

                    accessions_downloaded.append(accs)

    # Done!
    return accessions_downloaded


def get_cached_gene_seqs_by_id(gene_ids, download_missing=False):
    '''A generator method that provides the sequence of the requested genes.

    Arguments:
        gene_ids (list):         A list of gene IDs to use in the generator
        download_missing (bool): If a gene hasn't been cached, download it.

    '''

    if download_missing:
        to_download = []

        for gene_id in gene_ids:
            if not os.path.exists(get_gene_cache(gene_id)):
                to_download.append(gene_id)

        download_ncbi_genes(to_download)

    for gene_id in gene_ids:
        fna_fp = os.path.join(get_gene_cache(gene_id), f"{gene_id}.fna")

        with open(fna_fp, 'r') as fp:
            for header, seq in parse_fna(fp):
                yield gene_id, seq


def get_properties_from_ncbi_fasta_header(header, key=None):
    '''The NCBI Datasets API provides a multi-FASTA file when downloading genes.
    The headers are formatted like this:
        >NC_045512.2:26245-26472 E [organism=coronavirus 2] [gene_id=43740570] [chromosome=]
    This method will extract the value of the specified key from a FASTA header.

    In addition to the properties provided in square brackets, the following are also supported:
        - `accession`
        - `start`
        - `end`
        - `name`

    Arguments:
        header (string): The NCBI gene FASTA-header to parse
        key (string):    The key of the field to extract from the header.
                         If None, a dictionary of all properties is returned.

    Returns:
        See argument `key`.
    '''
    if header[0] == '>':
        header = header[1:]

    props = {}
    for i, prop in enumerate(header.split(' ')):
        if i == 0:
            # e.g., `NC_045512.2:26245-26472`
            acces, position = prop.split(':')
            start, end = position.split('-')
            props['accession'] = acces
            props['start'] = start
            props['end'] = end

        elif i == 1:
            # `E`
            props['name'] = prop

        else:
            # e.g., `[gene_id=43740570]`
            if prop[0] == '[' and prop[-1] == ']':
                k, v = prop[1:-1].split('=')
                props[k] = v

    if key is None:
        return props
    else:
        if key in props:
            return props[key]
        raise ValueError(f'Could not find `{key}` in header: `{header}`')


def get_assembly_cache(accession, mkdir=True):
    '''Generates a directory to store data for a particular accession.
    This function should be changed to your needs but keep in mind that changes
    may break the logic structure of your filesystem.

    Arguments:
        accession (string): The NCBI accession used for creating the directory
        mkdir (bool):       True, the directory is created.
                            False, the directory is not created.
                            Set to False if you just want to retrieve the path to the cache.

    Returns:
        The path to the directory where the NCBI data will be cached
    '''
    path = os.path.join(config.CACHE, accession[:config.CACHE_PREFIX_LENGTH], accession)
    if mkdir:
        os.makedirs(path, exist_ok=True)
    return path


def get_gene_cache(gene_id, mkdir=True):
    '''Generates a directory to store data for a particular gene.
    This function should be changed to your needs but keep in mind that changes
    may break the logic structure of your filesystem.

    Arguments:
        gene_id (string): The NCBI gene ID used for creating the directory
        mkdir (bool):     True, the directory is created.
                          False, the directory is not created.
                          Set to False if you just want to retrieve the path to the cache.

    Returns:
        The path to the directory where the NCBI data will be cached
    '''
    path = os.path.join(config.CACHE, f"gene-{gene_id}")
    if mkdir:
        os.makedirs(path, exist_ok=True)
    return path


def create_issl_indexes(accessions,
    bin_extract=None,
    bin_issl=None,
    stdout_extract=subprocess.DEVNULL,
    stdout_issl=subprocess.DEVNULL,
    stderr_extract=subprocess.DEVNULL,
    stderr_issl=subprocess.DEVNULL,
    force=False
):
    '''Extract offtargets and create ISSL index for each FNA file in the provided accession

    Arguments:
        accessions (list):     A list of accessions to process
        bin_extract (string):  A path to the extractOfftargets utility available in Crackling.
                               Set in configuration.
        bin_issl (string):     A path to the isslCreateIndex binary from Crackling.
                               Set in configuration.
        stdout_extract (file): File-like ojbect to redirect stdout for extractOfftargets to
        stdout_issl (file):    File-like ojbect to redirect stdout for isslCreateIndex to
        stderr_extract (file): File-like ojbect to redirect stderr for extractOfftargets to
        stderr_issl (file):    File-like ojbect to redirect stderr for isslCreateIndex to
        force (bool):          Rerun even if their output files already exist

    Returns:
        A list of accessions which an index was created for
    '''

    # Convert gene_ids to list if a string is provided
    if isinstance(accessions, str):
        accessions = [accessions]

    # Set defaults
    if bin_extract is None:
        bin_extract = config.BIN_EXTRACT
        
    if bin_issl is None:
        bin_issl = config.BIN_ISSL_IDX

    # remove accessions if the index already exists
    accessions = [
        accs
        for accs in accessions
        if len(
            glob(
                os.path.join(
                    get_assembly_cache(accs), 
                    '*.issl'
                )
            )
        ) == 0
    ]

    indexesCreateFor = []

    for accession in accessions:

        # Determine cache location
        cache = get_assembly_cache(accession)

        # Have all indexes been built successfully?
        success = True

        # Logic check: there should only be one *.fna per accession
        num_fna_files = len(glob(os.path.join(cache, '*.fna')))
        if num_fna_files != 1:
            raise RuntimeError((
                f'Found {num_fna_files} FNA files for {accession}. '
                f'There should be one.'
            ))

        # For the .fna file associated with the accession
        fna_path = glob(os.path.join(cache, '*.fna'))[0]

        fna = os.path.basename(fna_path)

        ots_file = os.path.join(cache, f"{fna}.offtargets.txt")
        issl_file = os.path.join(cache, f"{fna}.offtargets.issl")

        if not os.path.exists(ots_file) or force:

            # Extract off-targets
            cmd = [
                bin_extract,
                ots_file,
                fna_path
            ]

            try:
                subprocess.run(cmd,
                    stdout=stdout_extract,
                    stderr=stderr_extract
                )
            except Exception as e:
                if config.VERBOSE:
                    print(f'Failed to run: {" ".join(cmd)}')
                success &= False

        if not os.path.exists(issl_file) or force:
            # Create ISSL index
            cmd = [
                bin_issl,
                ots_file,
                '20',
                '8',
                issl_file
            ]
            try:
                subprocess.run(cmd,
                    stdout=stdout_issl,
                    stderr=stderr_issl
                )
            except Exception as e:
                if config.VERBOSE:
                    print(f'Failed to run: {" ".join(cmd)}')
                success &= False
        
        if success:
            indexesCreateFor.append(accession)
            
    return indexesCreateFor


def get_accessions_from_ncbi_table_export(
    filePath=config.NCBI_HUMAN_VIRUSES_TABLE, 
    **kwargs
):
    '''The 'NCBI Genome Information by Organism' table is a good way to obtain a filtered list of
    NCBI accessions. This method parses the CSV file that is provided when clicking 'Download' on
    this page https://www.ncbi.nlm.nih.gov/genome/browse/. Please ensure the 'Assembly' column is
    present before exporting.
    
    Arguments:
        filePath (string):  (optional) The path to the CSV file provided when Download is clicked on 
                            the webpage. By default, the precompiled list provided with ViralCut 
                            will be used.
        **kwargs:           Keyword arguments forwarded to `csv.reader()`
        
    Returns:
        A list of NCBI accessions found in the file
    
    '''
    
    accessions = []
    
    # Find the `Assembly` column then extract all values in it
    with open(filePath, 'r') as fp:
        fp_rdr = csv.reader(fp, **kwargs)
        idx_assembly = None
        
        for i, row in enumerate(fp_rdr):
            if i == 0:
                if 'Assembly' in row:
                    idx_assembly = row.index('Assembly')
                    continue
                else:
                    raise ValueError('The column `Assembly` must be present in the file, yet it could not be found.')
            
            accessions.append(row[idx_assembly])
    
    return accessions

def get_cached_gene_information_by_id(gene_id):
    '''Given some gene identifier, collect as properties about the gene. 
    Raises an exception if the gene does not exist in the cache.
    
    Arguments:
        gene_id (int):  The NCBI gene identifier
        
    Returns:
        A dictionary of properties
    '''
    dir = get_gene_cache(gene_id, mkdir=False)
    report_path = os.path.join(dir, 'data_report.json')
    if os.path.exists(report_path):
        with open(report_path, 'r') as fp:
            return json.loads(fp.readline())
    
    if config.VERBOSE:
        print(f'Could not find data report for gene #{gene_id}')
        
    return {}

def get_cached_tax_ids():
    '''This function crawls the assemblies cache to extract taxonomy IDs
    
    Returns:
        A list of taxonomy IDs
    '''
    
    tax_ids = []
    
    for path, subdirs, files in os.walk(config.CACHE):
        for name in files:
            filename = os.path.basename(name)
            if filename == 'data_report.json':
                with open(os.path.join(path, filename), 'r') as fp:
                    report = json.loads(fp.readline())
                if 'taxId' in report:
                    tax_ids.append(report['taxId'])

    return tax_ids
