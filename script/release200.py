# Dependencies
from dask.distributed import LocalCluster
from dask_jobqueue import LSFCluster
import argparse
import json
import sys
import os
import re

# Set path for custom dependencies
sys.path.append(os.path.dirname(__file__) + '/..')

# Custom dependencies
from src.pipeline.build import Build
from src.pipeline.build import get_cluster_names

# Constants
ROOT_PATH = os.path.dirname(__file__) + '/..'
LINCLUST_PATH = ROOT_PATH + '/data/clusters/chunk*.tsv.gz'
MGNIFAM_PATH = ROOT_PATH + '/data/mgnify/chunk*.fa.gz'
UNIPROT_PATH = ROOT_PATH + '/data/uniprot/chunk*.fa.gz'

# Main
if __name__ == '__main__':

    # Define new argument parser
    parser = argparse.ArgumentParser(description='Make new MGnifam release')
    # Feed command line arguments to it
    parser.add_argument(
        '-o', '--release_path', type=str, required=True,
        help='Path where to store new release'
    )
    parser.add_argument(
        '-i', '--input_path', nargs='+', type=str, required=True,
        help='Path to files holding cluster names'
    )
    parser.add_argument(
        '--shuffle', type=int, default=0,
        help='Whether to shuffle input cluster names'
    )
    parser.add_argument(
        '--batch_size', type=int, default=100,
        help='Maximum number of clusters to elaborate at each iterations'
    )
    parser.add_argument(
        '--max_clusters', type=int, default=0,
        help='Maximum number of clusters to elaborate'
    )
    parser.add_argument(
        '--linclust_path', type=str, default=LINCLUST_PATH,
        help='Path to LinClust clusters dataset'
    )
    parser.add_argument(
        '--uniprot_path', type=str, default=UNIPROT_PATH,
        help='Path to UniProt clusters dataset'
    )
    parser.add_argument(
        '--uniprot_height', type=int, required=False,
        help='Number of sequences in UniProt dataset'
    )
    parser.add_argument(
        '--uniprot_width', type=int, required=False,
        help='Length of longest sequence in UniProt dataset'
    )
    parser.add_argument(
        '--mgnifam_path', type=str, default=MGNIFAM_PATH,
        help='Path to MGnifam clusters dataset'
    )
    parser.add_argument(
        '--mgnifam_height', type=int, required=False,
        help='Number of sequences in MGnifam dataset'
    )
    parser.add_argument(
        '--mgnifam_width', type=int, required=False,
        help='Length of longest sequence in MGnifam dataset'
    )
    parser.add_argument(
        '--mobidb_cmd', nargs='+', type=str, default=['mobidb_lite.py'],
        help='Command for MobiDB Lite disorder predictor script'
    )
    parser.add_argument(
        '--muscle_cmd', nargs='+', type=str, default=['muscle'],
        help='Command for Muscle Multiple Sequence Alignment script'
    )
    parser.add_argument(
        '--hmmsearch_cmd', nargs='+', type=str, default=['hmmsearch'],
        help='Command for HMM search script'
    )
    parser.add_argument(
        '--hmmbuild_cmd', nargs='+', type=str, default=['hmmbuild'],
        help='Command for HMM build script'
    )
    parser.add_argument(
        '--hmmalign_cmd', nargs='+', type=str, default=['hmmalign'],
        help='Command for HMM align script'
    )
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Whether to print verbose log or not'
    )
    parser.add_argument(
        '-e', '--e_value', type=float, default=0.01,
        help='E-value threhsold for both UniProt and MGnifam comparisons'
    )
    parser.add_argument(
        '-c', '--cores', type=int, default=2,
        help='Maximum number of jobs per core'
    )
    parser.add_argument(
        '-m', '--memory', type=str, default='16 GB',
        help='Maximum amount of memory per job'
    )
    parser.add_argument(
        '-q', '--queue', type=str, required=None,
        help='Queue where to run jobs'
    )
    parser.add_argument(
        '--min_jobs', type=int, default=0,
        help='Minimum number of jobs to run in parallel'
    )
    parser.add_argument(
        '--max_jobs', type=int, default=100,
        help='Maximum number of jobs to run in parallel'
    )
    parser.add_argument(
        '--cluster_type', type=str, default='LSF',
        help='Type of cluster to use (only LSFCluster implemented yet)'
    )
    parser.add_argument(
        '--walltime', type=str, default='48:00',
        help='Time a job must be kept active, HH:MM format'
    )
    parser.add_argument(
        '--environ_path', type=str, required=False,
        help='Path to JSON storing environmental variables'
    )
    # Parse the arguments defined above
    args = parser.parse_args()

    # Define list of available cluster types
    cluster_types = {
        'LSF': LSFCluster,
        'local': LocalCluster
    }

    # Retrieve cluster type
    cluster_type = args.cluster_type
    # Set type of cluster
    if cluster_type not in cluster_types.keys():
        # Raise new exception and interrupt execution
        raise KeyError(' '.join([
            'Given cluster type {:s}'.format(cluster_type),
            'is not available, please choose one among these:',
            '[' + ', '.join([*cluster_types.values()]) + ']'
        ]))
    #
    # # Set environmental variables
    # env = {**os.environ.copy(), **{
    #     'PATH': ':'.join([
    #         "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/make",
    #         "/nfs/production/xfam/pfam/software/bin",
    #         "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/mgnifam",
    #         "/usr/lib64/qt-3.3/bin",
    #         "/ebi/lsf/ebi/ppm/10.2/bin",
    #         "/ebi/lsf/ebi/ppm/10.2/linux2.6-glibc2.3-x86_64/bin",
    #         "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/etc",
    #         "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/bin",
    #         "/usr/lpp/mmfs/bin",
    #         "/usr/local/bin",
    #         "/usr/bin",
    #         "/usr/local/sbin",
    #         "/usr/sbin",
    #         "/bin",
    #         "/usr/bin",
    #         "/homes/dclementel/bin"
    #     ]),
    #     'PYTHONPATH': ':'.join([
    #         '/ebi/sp/pro1/interpro/python-modules/lib64/python',
    #         '/ebi/sp/pro1/interpro/python-modules/lib/python'
    #     ])
    # }}

    # Initialize environmental variables
    env = dict()
    # Retrieve path to environment file
    environ_path = args.environ_path
    # Check if environment dictionary has been set
    if os.path.isfile(environ_path):
        # Open environmental variables file
        with open(environ_path, 'r') as env_file:
            # Update environmental variables
            env = json.load(env_file)
    # Loop through every key in retrieved dictionary
    for key, value in env.items():
        # Case calue is not a list
        if not isinstance(value, list):
            # Skip iteration, keep entry as is
            continue
        # Otherwise, join list as string, using double dots separators
        env[key] = ':'.join(value)
    # Merge current environmental variables with retrieved ones
    env = {**os.environ.copy(), **env}

    # Retrieve release path
    release_path = args.release_path
    # Check existence output directory
    if not os.path.exists(release_path):
        # Make output directory
        os.mkdir(release_path)

    # Initialize build pipeline
    pipeline = Build(
        # Pipeline parameters, required to handle job scheduling
        cluster_type=cluster_types.get(cluster_type),
        cluster_kwargs={
            'cores': args.cores,
            'memory': args.memory,
            'queue': args.queue,
            'walltime': args.walltime,
            'use_stdin': True,
            'processes': 1
        },
        # Path to datasets (fixed)
        linclust_path=args.linclust_path,
        mgnifam_path=args.mgnifam_path,
        uniprot_path=args.uniprot_path,
        # Compositional bias threshold settings
        comp_bias_threshold=0.2, comp_bias_inclusive=True,
        # Automatic trimming settings
        trim_threshold=0.4, trim_inclusive=True,
        filter_threshold=0.5, filter_inclusive=True,
        # Post trimming settings
        seed_min_width=1, seed_min_height=1,
        # Search against UniProt settings
        uniprot_e_value=args.e_value,
        uniprot_height=args.uniprot_height,
        uniprot_width=args.uniprot_width,
        # Search against MGnifam settings
        mgnifam_e_value=args.e_value,
        mgnifam_height=args.mgnifam_height,
        mgnifam_width=args.mgnifam_width,
        # Command line arguments
        mobidb_cmd=args.mobidb_cmd,  # Path to MobiDB Lite predictor
        muscle_cmd=args.muscle_cmd,  # Path to Muscle alignment algorithm
        hmm_build_cmd=args.hmmbuild_cmd,  # Path to hmmbuild script
        hmm_search_cmd=args.hmmsearch_cmd,  # Path to hmmsearch script
        hmm_align_cmd=args.hmmalign_cmd,  # Path to hmmalign script
        # Environmental variables
        env=env
    )

    # Retrieve input path and shuffle selector
    input_path = args.input_path
    shuffle = bool(args.shuffle)
    # Define clusters iterator
    cluster_names = get_cluster_names(paths=input_path, shuffle=shuffle)

    # Define number of clusters
    num_clusters = len(cluster_names)
    # Define maximum number of clusters to elaborate
    max_clusters = args.max_clusters if args.max_clusters else num_clusters
    # Threshold maximum number of clusters to elaborate
    num_clusters = min(num_clusters, max_clusters)
    # Retrieve batch size
    batch_size = args.batch_size
    # Loop through each batch of cluster names
    for i in range(0, num_clusters, batch_size):
        # Define batch index
        batch_index = str(i // batch_size)
        # Define directory for current batch
        batch_path = os.path.join(release_path, 'batches', batch_index)
        # make directory for current batch
        os.makedirs(batch_path, exist_ok=True)
        # Run pipeline for current batch of cluster names
        pipeline(
            cluster_names=cluster_names[i:min(num_clusters, i+batch_size)],
            clusters_path=batch_path,
            min_jobs=args.min_jobs,
            max_jobs=args.max_jobs,
            verbose=bool(args.verbose)
        )
