# Dependencies
from dask_jobqueue import LSFCluster
from glob import glob
import argparse
import sys
import os

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')

# Custom dependencies
from src.pipeline.hmm import HMMPipeline


# # Define root foldepath
# ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# # Define output path
# OUT_PATH = ROOT_PATH + '/tmp/seed'
# # Define path to clusters file
# CLUSTERS_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
# # Define path to MGnify sequences
# MGNIFY_PATH = ROOT_PATH + '/data/mgnify/chunk*.fa.gz'

# Main
if __name__ == '__main__':

    # Read command line arguments
    parser = argparse.ArgumentParser(
        description='Prepare dataset for releasing MGnify'
    )
    parser.add_argument(
        '--clusters_dir', type=str, default='./',
        help='Path to directory where MGYP* clusters are stored'
    )
    parser.add_argument(
        '--target_path', nargs='+', type=str, required=True,
        help='Path to target dataset (can be either .fasta or .fasta.gz)'
    )
    parser.add_argument(
        '-e', '--e_value', type=float, default=0.01,
        help='E-value threshold, only matches above will be considered'
    )
    parser.add_argument(
        '-z', '--z_score', type=float, required=False,
        help='Z-score: total length of the dataset, computed if not set'
    )
    parser.add_argument(
        '-c', '--cores', type=int, default=1,
        help='Number of cores per job'
    )
    parser.add_argument(
        '-m', '--memory', type=str, default='2 GB',
        help='Memory allocated per job'
    )
    parser.add_argument(
        '--walltime', type=str, default='48:00',
        help='Maximum time a job must be kept alive'
    )
    parser.add_argument(
        '--min_jobs', type=int, default=1,
        help='Minimum number of jobs to keep alive'
    )
    parser.add_argument(
        '--max_jobs', type=int, default=100,
        help='Maximum number of jobs to keep alive'
    )
    parser.add_argument(
        '--log_path', type=str, default='',
        help='Path where to store log (JSON)'
    )
    args = parser.parse_args()

    # Initialize cluster names
    cluster_names = list()
    # Loop through each cluster name in given folder
    for cluster_path in glob(os.path.join(args.clusters_dir, 'MGYP*')):
        # Get cluster name
        cluster_name = os.path.basename(cluster_path)
        # Store cluster name
        cluster_names.append(cluster_name)

    # Initialize pipeline
    pipeline = HMMPipeline(
        # Set type of cluster
        cluster_type=LSFCluster,
        # Set arguments for instantiate new cluster
        cluster_kwargs={
            'cores': args.cores,  # Maximum number of cores per job
            'processes': 1,  # Maximum number of processes
            'memory': args.memory,  # Memory allocated per job
            'walltime': args.walltime,  # Time before shutting down worker jobs
            'use_stdin': True  # Pass commands as stdin
        }
    )

    # # Define path to output directory
    # clusters_dir=OUT_PATH
    # Run pipeline
    pipeline(
        # Set cluster names
        cluster_names=cluster_names,
        # Set clusters directory
        clusters_dir=args.clusters_dir,
        # Set target dataset path
        target_path=args.target_path,
        # Set e-value threshold
        e_value=args.e_value,
        # Set z_score
        z_score=args.z_score,
        # Set minimum/maximum number of jobs
        min_jobs=args.min_jobs,
        max_jobs=args.max_jobs,
        # Path where to store log
        log_path=args.log_path
    )
