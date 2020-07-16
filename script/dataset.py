# Dependencies
from dask_jobqueue import LSFCluster
from dask.distributed import Client
import argparse
import time
import sys
import os

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__) + '/..'))

# Custom dependencies
from src.pipeline.dataset import DatasetPipeline
from src.dataset import LinClust, Fasta


# Main
if __name__ == '__main__':

    # Read command line arguments
    parser = argparse.ArgumentParser(
        description='Prepare dataset for releasing MGnify'
    )
    parser.add_argument(
        '-lp', '--linclust_path', type=str, default='',
        help='Path to LinClust results dataset (TSV)'
    )
    parser.add_argument(
        '--linclust_out_path', type=str, default='data/clusters',
        help='Path to LinClust chunks, must contain {:d} to enable indexing'
    )
    parser.add_argument(
        '--linclust_chunk_size', type=float, default=1e07,
        help='Maximum size of LinClust chunks'
    )
    parser.add_argument(
        '-up', '--uniprot_path', type=str, default='',
        help='Path to UniProt dataset (FASTA)'
    )
    parser.add_argument(
        '--uniprot_out_path', type=str, default='data/uniprot',
        help='Path to UniProt chunks, must contain {:d} to enable indexing'
    )
    parser.add_argument(
        '--uniprot_chunk_size', type=float, default=1e07,
        help='Maximum size of UniProt chunks'
    )
    parser.add_argument(
        '-mp', '--mgnify_path', type=str, default='',
        help='Path to MGnify dataset (FASTA)'
    )
    parser.add_argument(
        '--mgnify_out_path', type=str, default='data/mgnify',
        help='Path to UniProt chunks, must contain {:d} to enable indexing'
    )
    parser.add_argument(
        '--mgnify_chunk_size', type=float, default=1e07,
        help='Maximum size of MGnify chunks'
    )
    parser.add_argument(
        '-c', '--cores', type=int, default=1,
        help='Number of cores per job'
    )
    parser.add_argument(
        '-m', '--memory', type=str, default='16 GB',
        help='Memory allocated per job'
    )
    parser.add_argument(
        '--walltime', type=str, default='48:00',
        help='Maximum time a job must be kept alive'
    )
    parser.add_argument(
        '--log_path', type=str, default='',
        help='Path where to store log (JSON)'
    )
    args = parser.parse_args()


# # Define script root path
# ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# # Define path to clusters file
# CLUSTERS_PATH = '/nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz'
# # Define path to UniProt sequneces
# UNIPROT_PATH = '/nfs/production/xfam/pfam/data/pfamseq/uniprot'
# # Define path to MGnify sequences
# MGNIFY_PATH = '/nfs/production/xfam/pfam/data/mgnify/mgnify.fa.gz'
# # Define dataset path
# DATA_PATH = ROOT_PATH + '/data'


# Initialize dataset pipeline
pipeline = DatasetPipeline(
    # Define cluster type
    cluster_type=LSFCluster,
    # Define cluster parameters
    cluster_kwargs={
        'cores': args.cores,
        'memory': args.memory,
        'walltime': args.walltime,
        'use_stdin': True
    }
)

# Define dataset paths
linclust_out_dir = os.path.dirname(args.linclust_out_path)
uniprot_out_dir = os.path.dirname(args.uniprot_out_path)
mgnify_out_dir = os.path.dirname(args.mgnify_out_path)

# Initialize LinClust directory
if not os.path.exists(linclust_out_dir):
    # Try to make LinClust directory
    os.mkdir(linclust_out_dir)

# Initialize UniProt directory
if not os.path.exists(uniprot_out_dir):
    # Try to make UniProt directory
    os.mkdir(uniprot_out_dir)

# Initialize MGnify directory
if not os.path.exists(mgnify_out_dir):
    # Try to make MGnify directory
    os.mkdir(mgnify_out_dir)

# Run the pipeline
pipeline(
    # LinClust dataset
    linclust_path=args.linclust_path,
    linclust_chunk_size=args.linclust_chunk_size,
    linclust_out_path=args.linclust_out_path,
    # UniProt dataset
    uniprot_path=args.uniprot_path,
    uniprot_chunk_size=args.uniprot_chunk_size,
    uniprot_out_path=args.uniprot_out_path,
    # MGnify path
    mgnify_path=args.mgnify_path,
    mgnify_chunk_size=args.mgnify_chunk_size,
    mgnify_out_path=args.mgnify_out_path,
    # Define log path
    log_path=args.log_path
)

# Run the pipeline
# # Initialize new cluster
# cluster = LSFCluster(
#     cores=1,  # Number of cores per job
#     memory='8GB',  # Memory allocated per job
#     walltime='48:00',  # Time before shutting down worker jobs
#     log_directory=OUT_PATH+'/tmp',  # Logs directory
#     use_stdin=True  # Pass commands as stdin
# )
# # Require two jobs, one for Clusters, one for MGnify
# cluster.scale(jobs=2)
# # Debug
# print(cluster)
#
# # Initialize new client
# client = Client(cluster)
# # Debug
# print(client)
#
# # Debug
# print('Chunking clusters dataset...')
# # Make new instance of clusters dataset
# ds_cluster = Cluster(path=CLUSTERS_PATH)
# # Submit chunking function using bsub
# ft_cluster = client.submit(
#     # Chunking function to submit
#     ds_cluster.to_chunks,
#     # Formattable path of each chunk
#     chunk_path=OUT_PATH+'/data/clusters/chunk{:04d}.tsv.gz',
#     # Size (number of lines) of every cluster
#     chunk_size=1e07
# )
#
# # Debug
# print('Chunking MGnify dataset...')
# # Make new instance of MGnify dataset
# ds_mgnify = MGnify(path=MGNIFY_PATH)
# # Submit chunking function using bsub
# ft_mgnify = client.submit(
#     # Chunking function to submit
#     ds_mgnify.to_chunks,
#     # Formattable path of each chunk
#     chunk_path=OUT_PATH+'/data/mgnify/chunk{:04d}.fa.gz',
#     # Size (number of lines) of every cluster
#     chunk_size=1e06
# )
#
# # Initialize timers
# time_beg, time_end = time.time()
#
# # Gather futures from distributed memory
# client.gather([ft_cluster, ft_mgnify])
#
# # Release workers
# cluster.scale(jobs=0)
#
# # Update timers
# time_end = time.time()
# time_took = time_end - time_beg
# # Debug
# print('Took {:.0f} seconds to chunk both datasets'.format(time_took))
