# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client
import time
import sys
import os

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__) + '/..'))

# Custom dependencies
from src.dataset import Cluster
from src.dataset import MGnify


# Define path to clusters file
CLUSTERS_PATH = '/nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = '/nfs/production/xfam/pfam/data/mgnify/mgnify.fa.gz'
# Define script output path
OUT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'

# Initialize new cluster
cluster = LSFCluster(
    cores=1,  # Number of cores per job
    memory='8GB',  # Memory allocated per job
    walltime='48:00',  # Time before shutting down worker jobs
    log_directory=OUT_PATH+'/tmp',  # Logs directory
    use_stdin=True  # Pass commands as stdin
)
# Require two jobs, one for Clusters, one for MGnify
cluster.scale(jobs=2)
# Debug
print(cluster)

# Initialize new client
client = Client(cluster)
# Debug
print(client)

# Debug
print('Chunking clusters dataset...')
# Make new instance of clusters dataset
ds_cluster = Cluster(path=CLUSTERS_PATH)
# Submit chunking function using bsub
ft_cluster = client.submit(
    # Chunking function to submit
    ds_cluster.to_chunks,
    # Formattable path of each chunk
    chunk_path=OUT_PATH+'/data/clusters/chunk{:04d}.tsv.gz',
    # Size (number of lines) of every cluster
    chunk_size=1e07
)

# Debug
print('Chunking MGnify dataset...')
# Make new instance of MGnify dataset
ds_mgnify = MGnify(path=MGNIFY_PATH)
# Submit chunking function using bsub
ft_mgnify = client.submit(
    # Chunking function to submit
    ds_mgnify.to_chunks,
    # Formattable path of each chunk
    chunk_path=OUT_PATH+'/data/mgnify/chunk{:04d}.fa.gz',
    # Size (number of lines) of every cluster
    chunk_size=1e06
)

# Initialize timers
time_beg = time.time()
time_end = None
time_took = None

# Gather futures from distributed memory
client.gather([ft_cluster, ft_mgnify])

# Release workers
cluster.scale(jobs=0)

# Update timers
time_end = time.time()
time_took = time_end - time_beg
# Debug
print('Took {:.0f} seconds to chunk both datasets'.format(time_took))
