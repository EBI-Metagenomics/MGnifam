# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client
import time

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

# Initialize timers
time_beg = time.time()
time_end = None
time_took = None

# Debug
print('Chunking {:s} with chunk size {:d}'.format(CLUSTERS_PATH, int(1e07)))
# Make chunking
future = client.submit(
    ds.chunk_clusters,
    out_path=OUT_PATH + '/data/clusters/{:06d}.tsv.gz',
    chunk_size=int(1e07)
)
# Get the results
_ = future.result()

# Update timers
time_end = time.time()
time_took = time_end - time_beg
# Debug
print('Took {:.0f} seconds to chunk {:s}'.format(time_took, CLUSTERS_PATH))

# # Chunk mgnify file
# ds.chunk_mgnify(
#     out_path=OUT_PATH + '/data/mgnify/{}.fs.gz',
#     chunk_size=int(1e06)
# )
