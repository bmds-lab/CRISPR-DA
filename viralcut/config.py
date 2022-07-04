# The root taxonomy ID used. By default, the viruses superkingdom (ID 10239)
# is the root.
ROOT_TAX_ID = 10239

# The directory to write files downloaded from NCBI to
CACHE = '/mnt/ssd1/genomes-viruses'

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
BIN_ISSL_SCORE = '/home/jake/Crackling-v2.0/bin/isslScoreOfftargets'
BIN_ISSL_SCORE = '/mnt/hdd1/CRISPR/genomes/S-aureus-assemblies/isslScoreOfftargets_8b36d93'

# Path to the ISSL indexing binary
BIN_ISSL_IDX = '/home/jake/Crackling-v2.0/bin/isslCreateIndex'

# Path to the CRISPR site extraction utility available in Crackling
BIN_EXTRACT = 'extractOfftargets'

