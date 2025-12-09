from importlib import resources

# The 'NCBI Genome Information by Organism' table is a good way to obtain a filtered list of
# NCBI accessions. This method parses the CSV file that is provided when clicking 'Download' on
# this page https://www.ncbi.nlm.nih.gov/genome/browse/. Please ensure the 'Assembly' column is
# present before exporting. 
NCBI_HUMAN_VIRUSES_TABLE = resources.path('crispr_da.resources', 'human_viruses.csv')

# The directory to write files downloaded from NCBI to
CACHE = './genomes-viruses'

# When accessing the NCBI Datasets via the REST API, the requested URL may 
# rejected due to exceeding a reasonable length. To avoid this, queires 
# will be batched by this size.
NCBI_BATCH_SIZE = 100

# Print messages during processing
VERBOSE = True

# The consensus value used in Crackling
CONSENSUS_N = 2

# Path to the ISSL scoring binary
BIN_ISSL_SCORE = 'isslScoreOfftargets'

# Path to the ISSL indexing binary
BIN_ISSL_IDX = 'isslCreateIndex'

# Path to the CRISPR site extraction utility available in Crackling
BIN_EXTRACT = 'extractOfftargets'