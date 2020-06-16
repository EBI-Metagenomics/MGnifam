# Common dependencies
import numpy as np
import subprocess
import argparse
import tempfile
import time
import sys
import os
import re

# Updating python path
sys.path.append(os.path.dirname(os.path.realpath(__file__) + '/..'))

# Custom dependencies
from src.msa import MSA
from src.bjob import Bjob
# from ..src.transform import Compose
# from ..src.transform import OccupancyTrim
# from ..src.transform import OccupancyFilter
# from ..src.transform import MakeNonRedundant


# Read arguments
parser = argparse.ArgumentParser(description='Develop a new MGnifam release')
parser.add_argument(
    '-i', '--in_clusters', type=str, nargs='+',
    help='Name of input clusters, returned by LinClust'
)
parser.add_argument(
    '-o', '--out_dir', type=str, required=True,
    help='Path to directory where outcomes will be stored'
)
parser.add_argument(
    '--batch_size', type=int, default=1000,
    help='Number of clusters which can be processed at the same time'
)
args = parser.parse_args()

# Get output directory path
out_dir = re.sub(r'/$', '', args.out_dir)

# Get input clusters names
in_clusters = args.in_clusters
# Get number of clusters
num_clusters = len(in_clusters)
# Get batch size
batch_size = args.batch_size
# # Define index of last batch
# last_batch = (num_clusters // batch_size)
# # Define number of batches
# num_batches = 1 + last_batch

# Run mgseed.pl for each batch of clusters
for i in range(0, num_clusters, batch_size):
    # Get current batch of cluster names
    batch_clusters = in_clusters[i:min(i+batch_size, num_clusters)]
    # Define batch path (make temporary directory)
    batch_path = tempfile.mkdtemp()
    # Define set running jobs
    running = set()
    # Debug
    print('Running mgseed.pl for all the {} clusters in current batch'.format(
        len(batch_clusters)
    ))

    # Loop through each cluster in current batch
    for cluster_name in batch_clusters:
        # Run mgseed.pl in current batch directory
        out = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            cwd=batch_path,  # Set directory
            args=['mgseed.pl', '-cluster', cluster_name]
        )
        # Debug
        print('mgseed.pl:', out)
        # Get process id as string
        job_id = Bjob.id_from_string(out.stdout)
        # Debug
        print('Retrieved job id: {}'.format(job_id))
        # Save job id
        running.add(Bjob(id=job_id, status='RUN'))

    # Loop through each job while it is still running
    while running:
        # Get list of running jobs
        bjobs = list(running)
        # Query for job statuses
        is_running = list(map(lambda bjob: bjob.is_running(), bjobs))
        # Update running jobs
        running = set([bjobs[i] for i in range(len(is_running)) if is_running[i]])
        # Debug
        print('There are {} jobs which are still running:\n{}'.format(
            len(running),  # Number of running jobs
            ', '.join([job.id for job in running])  # Actual ids of running jobs
        ))
        # Set some delay (30 sec)
        time.sleep(30)

    # Run check_uniprot.pl in current batch directory
    out = subprocess.run(
        capture_output=True,
        encoding='utf-8',
        cwd=batch_path,
        args=[
            'check_uniprot.pl'
        ])
    # Debug
    print('check_uniprot.pl:', out)
