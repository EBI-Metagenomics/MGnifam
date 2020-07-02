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
import re

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')

# Custom dependencies
from src.pipeline import Seed


# Define root foldepath
ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# Define path to clusters file
CLUSTERS_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = ROOT_PATH + '/data/mgnify/chunk*.fa.gz'


if __name__ == '__main__':

    # Initialize distributed cluster
    cluster = LSFCluster(
        cores=1,  # Number of cores per job
        # scheduler_port=39156,  # Port where scheduler will listen
        memory='4GB',  # Memory allocated per job
        walltime='04:00',  # Time before shutting down worker jobs
        # log_directory=ROOT_PATH+'/tmp/log',  # Logs directory
        use_stdin=True  # Pass commands as stdin
    )
    # # Initialize local cluster
    # cluster = LocalCluster()
    cluster.adapt(minimum=20, maximum=100)
    # Debug
    print('Cluster:', cluster)
    print('Dashboard is at', cluster.dashboard_link)

    # Initialize new client
    client = Client(cluster)
    # Debug
    print('Client:', client)

    # Instantiate new seed alignment pipeline
    pipeline = Seed(
        clusters_path=CLUSTERS_PATH,
        mgnify_path=MGNIFY_PATH,
        dask_client=client
    )

    # # Define a single cluster name
    # cluster_names = ['MGYP001051398202']
    # # Define list of cluster names
    # cluster_names = [
    #     'MGYP000031769986', 'MGYP000002872901', 'MGYP000016484927',
    #     'MGYP000029146431', 'MGYP000029539082', 'MGYP000036242905',
    #     'MGYP000039062111', 'MGYP000040208386', 'MGYP000043559380',
    #     'MGYP000046805896', 'MGYP000055620076', 'MGYP000068575450',
    #     'MGYP000079685113', 'MGYP000089028905', 'MGYP000112885575',
    #     'MGYP000116866465', 'MGYP000124177638', 'MGYP000128352270',
    #     'MGYP000130617288'
    # ]
    # Initialize list of cluster names
    cluster_names = list()
    # Input file path
    in_path = '/hps/nobackup2/production/metagenomics/dclementel/MGnifam_build/Batch_lists/build_list0300'
    # Open input file
    with open(in_path, 'r') as in_file:
        # Read each file line
        for line in in_file:
            # Get cluster name
            match = re.search(r'^(\S+)', line)
            # Case line format does not match the expected one
            if not match:
                continue  # Skip iteration
            # Otherwise, add cluster name to input list
            cluster_names.append(match.group(1))
            # Early stopping
            if len(cluster_names) == 100:
                break

    # Define output path
    cluster_out_dir = ROOT_PATH + '/tmp'
    # Run the pipeline
    pipeline(
        cluster_names=cluster_names,
        cluster_dir=cluster_out_dir,
        verbose=True,
        log=True
    )
