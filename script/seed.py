"""Seed creation
Previously "create_seed_from_cluster.pl"
Given cluster names, retrieves alignment and computes compositional BIAS. If
compositional BIAS is above a certain threshold (e.g. 0.2), then discards the
current sequence alignment
"""

# Current cluster file is at:
# /nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz
# It is opened and unzipped every time there is the need of getting clusters sequences
# This can be greatly improved by leveraging distributed cluster features and
# Loading the whole dataset could be a good way of

# Dependencies
from dask.distributed import Client
# from dask.distributed import LocalCluster
from dask_jobqueue import LSFCluster
import time
import sys
import os

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')

# Custom dependencies
from src.pipeline import Seed


# Define root foldepath
ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# Define path to clusters file
CLUSTERS_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = ROOT_PATH + '/data/mgnify/chunk*.tsv.gz'


if __name__ == '__main__':

    # Initialize distributed cluster
    cluster = LSFCluster(
        cores=1,  # Number of cores per job
        memory='2GB',  # Memory allocated per job
        walltime='02:00',  # Time before shutting down worker jobs
        log_directory=ROOT_PATH+'/tmp',  # Logs directory
        use_stdin=True  # Pass commands as stdin
    )
    # # Initialize local cluster
    # cluster = LocalCluster()
    cluster.adapt(maximum_jobs=100)
    # Debug
    print('Cluster:', cluster)

    # Initialize new client
    client = Client(cluster)
    # Debug
    print('Client:', client)

    # Instantiate new seed alignment pipeline
    pipe = Seed(
        clusters_path=CLUSTERS_PATH,
        mgnify_path=MGNIFY_PATH,
        dask_client=client
    )

    # Define a single cluster name
    cluster_names = ['MGYP000848664103']
    # Run the pipeline
    pipe(cluster_names=cluster_names)
