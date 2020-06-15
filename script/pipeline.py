# # Custom dependencies
# from ..src.msa import MSA
# from ..src.transform import Compose
# from ..src.transform import OccupancyTrim
# from ..src.transform import OccupancyFilter
# from ..src.transform import MakeNonRedundant

# Common dependencies
import subprocess
import argparse
import tempfile
import os
import re


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
# Define number of batches
num_batches = 1 + (num_clusters // batch_size)
# Run Mgseed.pl for each batch of clusters
for i in range(0, num_batches):
    # Get current batch of cluster names
    batch_clusters = in_clusters[(i*batch_size):min((i+1)*batch_size, num_clusters)]
    # Define batch path (make temporary directory)
    batch_path = tempfile.mkdtemp()
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
        print(out)
    # TODO Run check_uniprot.pl in current batch directory
