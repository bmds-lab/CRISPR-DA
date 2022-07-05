from importlib import resources

# The 'NCBI Genome Information by Organism' table is a good way to obtain a filtered list of
# NCBI accessions. This method parses the CSV file that is provided when clicking 'Download' on
# this page https://www.ncbi.nlm.nih.gov/genome/browse/. Please ensure the 'Assembly' column is
# present before exporting. 
NCBI_HUMAN_VIRUSES_TABLE = resources.path('viralcut.resources', 'human_viruses.csv')


# The root taxonomy ID used. By default, the viruses superkingdom (ID 10239)
# is the root.
ROOT_TAX_ID = 10239

# The directory to write files downloaded from NCBI to
CACHE = '~./genomes-viruses'

# When downloading through NCBI Datasets PyLib, the requested URL may be
# rejected due to exceeding a reasonable length, if too many accessions are 
# requested at once. To avoid this, downloads will be batched by this size.
NCBI_ACCESSION_BATCH_SIZE = 100

# When downloading NCBI assemblies, the default cacheing method creates
# sub-directories named as the first few characters of the accession.
# Here, specify the number of characters to use.
# Example: `GCA_000841585.1` is cached in `GCA_000841/GCA_000841585.1`
CACHE_PREFIX_LENGTH = 10

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

