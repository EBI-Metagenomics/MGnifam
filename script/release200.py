# Dependencies
from dask.distributed import LocalCluster
from dask_jobqueue import LSFCluster
import argparse
import sys
import os
import re

# Set path for custom dependencies
sys.path.append(os.path.dirname(__file__) + '/..')

# Custom dependencies
# from src.pipeline.release import ReleasePipeline
from src.pipeline.batch import Batch
from src.utils import get_paths

# Constants
ROOT_PATH = os.path.dirname(__file__) + '/..'
LINCLUST_PATH = ROOT_PATH + '/data_/clusters/chunk*.tsv.gz'
MGNIFAM_PATH = ROOT_PATH + '/data_/mgnify/chunk*.fa.gz'
UNIPROT_PATH = ROOT_PATH + '/data_/uniprot/chunk*.fa.gz'

# Iterate clusters file
def iter_clusters(in_path, batch_size=100, max_clusters=None):
    """Iterate through clusters files

    Args
    clusters_paths (str/list)   Path or list of paths to linclust files
    batch_size (int)            Maximum number of clusters per batch
    max_clusters (int)          Maximum number of clusters to iterate

    Return
    (generator)                 A generator that yelds batch as
                                tuple(batch index, cluster names)
    """
    # Initialize current cluster index
    cluster_index = 0
    # Initialize current batch index
    batch_index = 0
    # Initialize current batch of clusters
    batch_clusters = list()
    # Loop through each file input path
    for curr_path in get_paths(in_path):
        # Open file at current path
        with open(curr_path, 'r') as file:
            # Loop through each line in current linclust file path
            for line in file:

                # Check if current line is matches expected format
                match = re.search(r'^(\S+)\s+(\d+)', line)
                # Case line does not match expected format: skip
                if not match:
                    continue

                # Otherwise, retrieve cluster name and cluster size
                cluster_name = str(match.group(1))
                cluster_size = int(match.group(2))
                # Add current cluster name to batch
                batch_clusters.append(cluster_name)
                # Update current cluster index
                cluster_index += 1
                # Define batch index
                batch_index = (cluster_index - 1) // batch_size

                # Case cluster index has reached maximum size, exit loop
                if max_clusters is not None:
                    if cluster_index >= max_clusters:
                        break
                # Ensure cluster index has not reached batch size
                if (cluster_index % batch_size) != 0:
                    continue

                # Otherwise, yield batch index and list of clusters
                yield batch_index, batch_clusters
                # Reset batch of clusters
                batch_clusters = list()

            # Case cluster index has reached maximum size, exit loop
            if max_clusters is not None:
                if cluster_index >= max_clusters:
                    break

    # In case we reached this point, check for non returned clusters
    if len(batch_clusters) > 0:
        # Yield last batch
        yield batch_index, batch_clusters

# Main
if __name__ == '__main__':

    # Define list of available cluster types
    cluster_types = {
        'LSF': LSFCluster,
        'local': LocalCluster
    }

    # Define new argument parser
    parser = argparse.ArgumentParser(description='Make new MGnifam release')
    # Feed command line arguments to it
    parser.add_argument(
        '-o', '--out_path', type=str, required=True,
        help='Output path where release will be stored'
    )
    parser.add_argument(
        '-i', '--in_path', type=str, required=True,
        help='Path where files holding cluster names are stored'
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
        '--mgnifam_path', type=str, default=MGNIFAM_PATH,
        help='Path to MGnifam clusters dataset'
    )
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Whether to print verbose log or not'
    )
    parser.add_argument(
        '-c', '--cores', type=int, default=1,
        help='Maximum number of jobs per core'
    )
    parser.add_argument(
        '-m', '--memory', type=str, default='8 GB',
        help='Maximum amount of memory per job'
    )
    parser.add_argument(
        '--cluster', type=str, default='LSF',
        help='Type of cluster to use (only LSFCluster implemented yet)'
    )
    parser.add_argument(
        '--env', type=str, required=False,
        help='Path to JSON storing environmental variables'
    )
    # Parse the arguments defined above
    args = parser.parse_args()

    # Set type of cluster
    if args.cluster not in cluster_types.keys():
        # Raise new exception and interrupt execution
        raise KeyError(' '.join([
            'Given cluster type {:s}'.format(args.cluster),
            'is not available, please choose one among these:',
            '[' + ', '.join(cluster_types.values()) + ']'
        ]))

    # Set environmental variables
    env = {**os.environ.copy(), **{
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

    # Initialize batch pipeline
    pipeline = Batch(
        # Pipeline parameters, required to handle job scheduling
        cluster_type=cluster_types.get(args.cluster),
        cluster_kwargs={
            'cores': args.cores,
            'memory': args.memory
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
        uniprot_e_value=0.01, uniprot_z_score=None,
        # Command line arguments
        mobidb_cmd=['mobidb_lite.py'],  # Path to MobiDB Lite predictor
        muscle_cmd=['muscle'],  # Path to Muscle alignment algorithm
        hmm_build_cmd=['hmmbuild'],  # Path to hmmbuild script
        hmm_search_cmd=['hmmsearch'],  # Path to hmmsearch script
        # Environmental variables
        env=env
    )

    # Get cluster names
    batch_index, cluster_names = next(iter_clusters(
        in_path=args.in_path,
        batch_size=100,
        max_clusters=100
    ))

    # Check existence output directory
    if not os.path.exists(args.out_path):
        # Make output directory
        os.mkdir(args.out_path)

    # Run the pipeline
    pipeline(
        cluster_names=cluster_names,
        clusters_path=args.out_path,
        min_jobs=1,
        max_jobs=100,
        verbose=bool(args.verbose)
    )
