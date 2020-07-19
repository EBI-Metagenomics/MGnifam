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
# from dask.distributed import LocalCluster
from dask_jobqueue import LSFCluster
import argparse
import time
import sys
import os
import re

# Update python path for custom dependencies
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/..')

# Custom dependencies
from src.pipeline.seed import SeedPipeline


# Define root foldepath
ROOT_PATH = '/nfs/production/metagenomics/mgnifams/dclementel/MGnifam'
# Define path to clusters file
LINCLUST_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
# Define path to MGnify sequences
MGNIFY_PATH = ROOT_PATH + '/data/mgnify/chunk*.fa.gz'


# # Main
# if __name__ == '__main__':
#
#     # Initialize distributed cluster
#     cluster = LSFCluster(
#         cores=1,  # Number of cores per job
#         memory='8GB',  # Memory allocated per job
#         walltime='04:00',  # Time before shutting down worker jobs
#         use_stdin=True  # Pass commands as stdin
#     )
#     # Set number of jobs boundaries
#     cluster.adapt(minimum=20, maximum=100)
#     # Debug
#     print('Cluster:', cluster)
#     print('Dashboard is at', cluster.dashboard_link)
#
#     # Initialize new client
#     client = Client(cluster)
#     # Debug
#     print('Client:', client)
#
#     # Instantiate new seed alignment pipeline
#     pipeline = Seed(
#         linclust_path=CLUSTERS_PATH,
#         mgnify_path=MGNIFY_PATH,
#         dask_client=client
#     )
#
#     # # Define a single cluster name
#     # cluster_names = ['MGYP001051398202']
#     # # Define list of cluster names
#     # cluster_names = [
#     #     'MGYP000031769986', 'MGYP000002872901', 'MGYP000016484927',
#     #     'MGYP000029146431', 'MGYP000029539082', 'MGYP000036242905',
#     #     'MGYP000039062111', 'MGYP000040208386', 'MGYP000043559380',
#     #     'MGYP000046805896', 'MGYP000055620076', 'MGYP000068575450',
#     #     'MGYP000079685113', 'MGYP000089028905', 'MGYP000112885575',
#     #     'MGYP000116866465', 'MGYP000124177638', 'MGYP000128352270',
#     #     'MGYP000130617288'
#     # ]
#     # Initialize list of cluster names
#     cluster_names = list()
#     # Set upper bound to clusters number
#     cluster_number = 1000
#     # Input file path
#     in_path = '/hps/nobackup2/production/metagenomics/dclementel/MGnifam_build/Batch_lists/build_list0300'
#     # Open input file
#     with open(in_path, 'r') as in_file:
#         # Read each file line
#         for line in in_file:
#             # Get cluster name
#             match = re.search(r'^(\S+)', line)
#             # Case line format does not match the expected one
#             if not match:
#                 continue  # Skip iteration
#             # Otherwise, add cluster name to input list
#             cluster_names.append(match.group(1))
#             # Case limit of clusters has been reached
#             if cluster_number == len(cluster_names):
#                 break  # Exit cycle
#
#     # Define output path
#     cluster_out_dir = ROOT_PATH + '/tmp/seed'
#     # Run the pipeline
#     pipeline(
#         cluster_names=cluster_names,
#         clusters_dir=cluster_out_dir,
#         comp_bias_threshold=0.2,
#         comp_bias_inclusive=True,
#         verbose=True
#     )

# Main
if __name__ == '__main__':

    # Get current environment
    env = os.environ.copy()
    # Setup local environment
    env = {**env, **{
        'PATH': ':'.join([
            "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/make",
            "/nfs/production/xfam/pfam/software/bin",
            "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/mgnifam",
            "/usr/lib64/qt-3.3/bin",
            "/ebi/lsf/ebi/ppm/10.2/bin",
            "/ebi/lsf/ebi/ppm/10.2/linux2.6-glibc2.3-x86_64/bin",
            "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/etc",
            "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/bin",
            "/usr/lpp/mmfs/bin",
            "/usr/local/bin",
            "/usr/bin",
            "/usr/local/sbin",
            "/usr/sbin",
            "/bin",
            "/usr/bin",
            "/homes/dclementel/bin"
        ]),
        'PYTHONPATH': ':'.join([
            '/ebi/sp/pro1/interpro/python-modules/lib64/python',
            '/ebi/sp/pro1/interpro/python-modules/lib/python'
        ])
    }}

    # Read command line arguments
    parser = argparse.ArgumentParser(
        description='Prepare dataset for releasing MGnify'
    )
    parser.add_argument(
        '--in_clusters', nargs='+', type=str, required=True,
        help='Path to one or more files containing cluster names'
    )
    parser.add_argument(
        '--release_path', type=str, required=True,
        help='Output directory path'
    )
    parser.add_argument(
        '--linclust_path', type=str, default=LINCLUST_PATH,
        help='Path to .tsv directory where clusters are stored'
    )
    parser.add_argument(
        '--mgnify_path', nargs='+', type=str, default=MGNIFY_PATH,
        help='Path to target dataset (can be either .fasta or .fasta.gz)'
    )
    parser.add_argument(
        '--comp_bias_threshold', type=float, default=0.2,
        help='Threshold maximum compositional bias'
    )
    parser.add_argument(
        '--comp_bias_inclusive', type=int, default=1,
        help='Wether to include or not compositional bias threshold value'
    )
    parser.add_argument(
        '--trim_occupancy', type=float, default=0.4,
        help='Minimum occupancy threshold used in automatic trimming'
    )
    parser.add_argument(
        '--trim_inclusive', type=int, default=1,
        help='Wether to include or not values equal to threshold when trimming'
    )
    parser.add_argument(
        '--filter_occupancy', type=float, default=0.5,
        help='Minimum occupancy threshold used in automatic filtering'
    )
    parser.add_argument(
        '--filter_inclusive', type=int, default=1,
        help='Wether to include or not values equal to threshold when filtering'
    )
    parser.add_argument(
        '-n', '--num_clusters', type=int, default=1000,
        help='Maximum number of clusters to run through the pipeline'
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
    # Get first batch in clusters
    _, cluster_names = next(SeedPipeline.iter_clusters(
        clusters_path=args.in_clusters,
        batch_size=args.num_clusters,
        max_clusters=args.num_clusters
    ))

    # Initialize pipeline
    pipeline = SeedPipeline(
        # Define path to LinCLust .tsv dataset
        linclust_path=args.linclust_path,
        # Define path to MGnify .fa[.gz] dataset
        mgnify_path=args.mgnify_path,
        # Set type of cluster
        cluster_type=LSFCluster,
        # Set arguments for instantiate new cluster
        cluster_kwargs={
            'cores': args.cores,  # Maximum number of cores per job
            'processes': 1,  # Maximum number of processes
            'memory': args.memory,  # Memory allocated per job
            'walltime': args.walltime,  # Time before shutting down worker jobs
            'use_stdin': True  # Pass commands as stdin
        },
        # Set environmental variables
        env=env
    )

    # Run pipeline
    pipeline(
        # TODO Define cluster names
        cluster_names=cluster_names,
        # Define output directory
        clusters_dir=args.release_path,
        # Set compositional bias threshold
        comp_bias_threshold=args.comp_bias_threshold,
        comp_bias_inclusive=args.comp_bias_inclusive,
        # Set minimum/maximum number of jobs
        min_jobs=args.min_jobs,
        max_jobs=args.max_jobs,
        # Path where to store log
        log_path=args.log_path
    )
