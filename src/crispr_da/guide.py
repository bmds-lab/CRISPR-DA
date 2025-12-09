import pandas as pd

CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNKNOWN = '?'
CODE_ERROR = '!'

class Guide:
    def __init__(self, seq):
        self.seq = seq

        # properties of a guide
        self.props = {}

        # scores
        self.assembly_scores = {}
        # defaultdict(lambda : {
        #     'score_name' : [],
        #     'score' : [],
        #     'unique_sites' : [],
        #     'total_sites' : []
        #     })

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

    def add_assembly_score(self, accession, mit, cfd, unique_sites, total_sites):
        self.assembly_scores[accession] = {
            'mit' : '',
            'cfd' : '',
            'unique_sites' : '',
            'total_sites' : ''
            }
        self.assembly_scores[accession]['mit'] = mit
        self.assembly_scores[accession]['cfd'] = cfd
        self.assembly_scores[accession]['unique_sites'] = unique_sites
        self.assembly_scores[accession]['total_sites'] = total_sites

    def assembly_scores_to_dataframe(self):
        return pd.DataFrame(self.assembly_scores.items())
