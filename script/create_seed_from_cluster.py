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


# # Dependencies
# import gzip


# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client
from dask.distributed import progress
# from dask.delayed import delayed
# import dask.dataframe as dd
# import pandas as pd
# from tqdm import tqdm
import gzip
import glob
import re
import os


# # Load clusters as pandas dataframe
# def read_tsv(in_path):
#     # Debug
#     print('Reading ', in_path)
#     # Make and return dataframe
#     return pd.read_csv(
#         in_path,  # Path to chunk
#         sep='\t',
#         usecols=[0, 1],
#         names=['cluster_name', 'sequence_acc'],
#         compression='infer'
#     )

# Get cluster entries from chunk
def get_sequences(in_name, in_path):
    # Dependencies
    import gzip
    import re
    # Get cluster entries
    sequences = list()
    # Unzip
    with gzip.open(in_path, 'rt') as in_file:
        # Loop through every line in file
        for line in in_file:
            # Search entries in text line
            found = re.search(r'^([a-zA-Z0-9]+)[ \t]+([a-zA-Z0-9]+)', line)
            # Case line does not match format
            if not found:
                # Skip iteration
                continue
            # Get cluster name and sequence accession
            cluster_name, sequence_acc = found.group(1), found.group(2)
            # Case current cluster name matches the input one
            if cluster_name == in_name:
                # Store entry
                sequences.append((cluster_name, sequence_acc))
    # Return entries found
    return sequences

# TODO Read cluster from input path
# Define path to clusters gzipped file (formatted as <cluster name>\t<seq acc>)
CLUSTER_PATH = '/nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz'

# Create new cluster instance
cluster = LSFCluster(cores=1, memory='2 GB', use_stdin=True)
# Create new client instance (interacts job scheduler)
client = Client(cluster)
# Require 12 jobs
cluster.scale(12)

# Debug
print('Cluster', cluster)
print('Client', client)

# Get input path
in_path = CLUSTER_PATH
# Check if it is gzipped
is_gzip = bool(re.search(r'\.gz$', in_path))
# Debug
print('Is gzipped? {}!'.format('Yes' if is_gzip else 'No'))
# Define default file handler
file_handler = lambda in_path: open(in_path, 'r')
# Define gzipped file handler
if is_gzip:
    file_handler = lambda in_path: gzip.open(in_path, 'rt')

# Define temporary directory where to store chunks
cluster_dir = './tmp/clusters'
# # Make cluster directory
# os.makedirs(cluster_dir, exist_ok=True)
# # Define number of rows per chunk
# batch_size = int(1e08)
# # Open input file
# with file_handler(in_path) as in_file:
#     # Intialize chunk
#     curr_chunk = None
#     # Intialize line counter
#     curr_line = 0
#     # Loop through each input file line
#     for line in tqdm(in_file):
#         # # Match values on each row
#         # found = re.search(r'^([a-zA-Z0-9]+)\t([a-zA-Z0-9]+)', line)
#         # # Skip iteration if values not matching
#         # if not found:
#         #     continue
#         # Check if current line is at the beginning of a new batch
#         if not (curr_line % batch_size):
#             # Eventually close previous chunk
#             if curr_chunk:
#                 curr_chunk.close()
#             # Create new chunk path
#             chunk_path = '{:s}/chunk{:d}.tsv.gz'.format(
#                 cluster_dir,  # Directory where chunk will be stored
#                 curr_line // batch_size  # Number of chunk stored
#             )
#             # Open new chunk
#             curr_chunk = gzip.open(chunk_path, 'wt')
#         # # Get values as cluster name and sequence accession
#         # cluster_name, sequence_acc = found.group(1), found.group(2)
#         # Append to list of sequences
#         curr_chunk.write(line)
#         # Update line counter
#         curr_line += 1
#     # Close last opened chunk
#     curr_chunk.close()


# Define a set of futures
futures = client.map(
    # Function that searches for cluster MGYP000853368667
    lambda in_path: get_sequences(in_path=in_path, in_name='MGYP000853368667'),
    # Chunks
    glob.glob(cluster_dir + '/chunk*.tsv.gz')
)
# Check progress
progress(futures)
# Wait for results
sequences = [future.result() for future in futures]
# Concatenate retrieved lists into a single list
sequences = [
    sequences[i][j]
    for i in range(len(sequences))
    for j in range(len(sequences[i]))
]
# Debug
print('Retrieved sequences', sequences)
# # Create a list of lazy functions ready to return a pandas.DataFrame
# dfs = [delayed(read_tsv)(path) for path in glob.glob(cluster_dir + '/chunk*.tsv.gz')]
# # Using delayed, assemble the pandas.DataFrames into a dask.DataFrame
# ddf = dd.from_delayed(dfs)
