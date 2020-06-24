# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client

# Custom dependencies
from src.pipeline.dataset import Dataset

# Define path to clusters file
CLUSTERS_PATH = '/nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = '/nfs/production/xfam/pfam/data/mgnify/mgnify.fa.gz'
# Define script output path
OUT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'

# Initialize new cluster
cluster = LSFCluster(cores=1, memory='4GB', use_stdin=True)
# Debug
print(cluster)

# Initialize new client
client = Client(cluster)
# Debug
print(client)

ds = Dataset(
    clusters_path=CLUSTERS_PATH,
    mgnify_path=MGNIFY_PATH
)

# Require one job
cluster.scale(1)
# Debug
print(cluster)

# Make chunking
future = client.submit(
    ds.chunk_clusters,
    out_path=OUT_PATH + '/data/clusters/{}.tsv.gz',
    chunk_size=int(1e06)
)
# Get the results
future.result()

# # Chunk mgnify file
# ds.chunk_mgnify(
#     out_path=OUT_PATH + '/data/mgnify/{}.fs.gz',
#     chunk_size=int(1e06)
# )
