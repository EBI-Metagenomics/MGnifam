# Common dependencies
import subprocess
import argparse
import tempfile
import shutil
import glob
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
# Define temporary build directory
build_path = tempfile.mkdtemp()
# # Define index of last batch
# last_batch = (num_clusters // batch_size)
# # Define number of batches
# num_batches = 1 + last_batch

# Run mgseed.pl for each batch of clusters
for i in range(0, num_clusters, batch_size):
    # Get current batch of cluster names
    batch_clusters = in_clusters[i:min(i+batch_size, num_clusters)]
    # Define batch path (make temporary directory)
    batch_path = out_dir + '/batch_{}'.format(i // batch_size)
    # Make batch directory
    os.mkdir(batch_path)

    # Initialize set of running jobs
    bjobs = list()
    # Debug
    print('Running mgseed.pl for all the {} clusters in batch {}'.format(
        len(batch_clusters),  # Number of clusters
        batch_path  # Current batch path
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
        # Save job id
        bjobs.append(Bjob(id=job_id, status='RUN'))

    # Check running mgseed.pl jobs
    Bjob.check(bjobs, delay=30)

    # Define set of running jobs
    bjobs = list()
    # Retrieve cluster paths
    cluster_paths = glob.glob(batch_path + '/MGYP*')
    # Debug
    print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
    # Run pfbuild for every cluster in batch
    for cluster_path in cluster_paths:
        # Run pfbuild current cluster directory
        out = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            cwd=cluster_path,  # Set directory
            args=['pfbuild', '-withpfmake', '-db', 'uniprot']
        )
        # Debug
        print('pfbuild:', out)
        # Get process id as string
        job_id = Bjob.id_from_string(out.stdout)
        # Save job id
        bjobs.append(Bjob(id=job_id, status='RUN'))

    # Check running pfbuild jobs
    Bjob.check(bjobs, delay=30)

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

    # Define kept clusters (MGnifam)
    kept_clusters = glob.glob(batch_path + '/MGYP*')
    # Define discarded clusters
    bias_clusters = glob.glob(batch_path + '/BIAS/MGYP*')
    uni_clusters = glob.glob(batch_path + '/Uniprot/MGYP*')
    # Loop through folders not discarded from this batch
    for cluster_path in kept_clusters:
        # Get cluster name
        cluster_name = os.path.basename(cluster_path)
        cluster_name = os.path.dirname(cluster_name)
        # Move cluster to build folder
        shutil.move(cluster_path, build_path + '/' + cluster_name)
        # Debug
        print('Moved {} to {}'.format(cluster_path, build_path))
    # Debug
    print('There are {:d} ({:.03f}%) possible MGnifam clusters: {}'.format(
        # Number of possible MGnifam
        len(kept_clusters),
        # Rate of possible MGnifam clusters
        len(kept_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible MGnifam
        ', '.join(os.path.basename(path) for path in kept_clusters)
    ))
    # Debug
    print('There are {:d} ({:.03f}%) clusters in BIAS: {}'.format(
        # Number of possible Pfam
        len(bias_clusters),
        # Rate of possible Pfam clusters
        len(bias_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible Pfam
        ', '.join(os.path.basename(path) for path in bias_clusters)
    ))
    # Debug
    print('There are {:d} ({:.03f}%) clusters in Uniprot: {}'.format(
        # Number of possible Pfam
        len(uni_clusters),
        # Rate of possible Pfam clusters
        len(uni_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible Pfam
        ', '.join(os.path.basename(path) for path in uni_clusters)
    ))
