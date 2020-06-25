"""Seed creation
Given a cluster names, retrieves alignment and computes compositional BIAS. If
compositional BIAS is above a certain threshold (e.g. 0.2), then discards the
current sequence alignment
"""

# Current cluster file is at:
# /nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz
# It is opened and unzipped every time there is the need of getting clusters sequences
# This can be greatly improved by leveraging distributed cluster features and
# Loading the whole dataset could be a good way of

# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client
import time

# Custom dependencies
from src.dataset import Cluster
from src.dataset import MGnify


# Define root foldepath
ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# Define path to clusters file
CLUSTERS_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = ROOT_PATH + '/data/mgnify/chunk*.tsv.gz'


# Initialize new cluster
cluster = LSFCluster(
    cores=2,  # Number of cores per job
    memory='2GB',  # Memory allocated per job
    walltime='02:00',  # Time before shutting down worker jobs
    log_directory=OUT_PATH+'/tmp',  # Logs directory
    use_stdin=True  # Pass commands as stdin
)
# Adapt cluster
cluster.adapt(maximum_jobs=100)
# Debug
print(cluster)

# Initialize new client
client = Client(cluster)
# Debug
print(client)

# Load (list of) clusters dataset from .tsv.gz chunks
ds_cluster = Cluster.from_str(CLUSTERS_PATH)
# Debug
print(ds_cluster)

# Load (list of) sequences dataset from .fa.gz chunks
ds_mgnify = MGnify.from_str(MGNIFY_PATH)
# Debug
print(ds_mgnify)

# Define input cluster name
cluster_names = ['MGYP000848664103']

# Initialize timers
time_beg = time.time()
time_end = None
time_took = None

# Search cluster members
futures = client.map(
    # Function to distribute
    lambda ds: ds.search(cluster_names),
    # List of cluster datasets
    ds_cluster
)
# retrieve found clusters
cluster_members = client.gather(futures)

# Update timers
time_end = time.time()
time_took = time_end - time_beg
# Debug
print('Took {:.0f} seconds to search for cluster members'.format(time_took))
print('got', cluster_members)
