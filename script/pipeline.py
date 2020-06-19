# Common dependencies
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import subprocess
import argparse
import shutil
import glob
import sys
import os
import re

# Set plotting non iteractive
matplotlib.use('Agg')

# Updating python path
sys.path.append(os.path.dirname(os.path.realpath(__file__) + '/..'))

# Custom dependencies
from src.msa import MSA
from src.bjob import Bjob
from src.transform import Compose
from src.transform import OccupancyTrim
from src.transform import OccupancyFilter
from src.transform import MakeNonRedundant


# Run mgseed
def run_mgseed(cluster_name, cwd='./', env=os.environ.copy()):
    # Run mgseed.pl in current batch directory
    out = subprocess.run(
        capture_output=True,  # Capture console output
        encoding='utf-8',  # Set output encoding
        env=env,
        cwd=cwd,  # Set directory
        args=['mgseed.pl', '-cluster', cluster_name]
    )
    # Debug
    print('mgseed.pl:', out)
    # Get process id as string
    job_id = Bjob.id_from_string(out.stdout)
    # Return new Bjob instance
    return Bjob(id=job_id, status='RUN')


# Run pfbuild
def run_pfbuild(cluster_path, env=os.environ.copy()):
    # Run pfbuild current cluster directory
    out = subprocess.run(
        capture_output=True,  # Capture console output
        encoding='utf-8',  # Set output encoding
        env=env,  # Set custom environment
        cwd=cluster_path,  # Set directory
        args=['pfbuild', '-withpfmake', '-db', 'uniprot']
    )
    # Debug
    print('pfbuild:', out)
    # Get process id as string
    job_id = Bjob.id_from_string(out.stdout)
    # Return new Bjob instance
    return Bjob(id=job_id, status='RUN')


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
# Define build directory
build_path = out_dir + '/build'
# Make build directory
os.mkdir(build_path)

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


# Run mgseed.pl for each batch of clusters
for i in range(0, num_clusters, batch_size):
    # Get current batch of cluster names
    batch_clusters = in_clusters[i:min(i+batch_size, num_clusters)]
    # Define batch path (make temporary directory)
    batch_path = out_dir + '/batch_{}'.format(i // batch_size)
    # Make batch directory
    os.mkdir(batch_path)

    # Debug
    print('Running mgseed.pl for all the {} clusters in batch {}'.format(
        len(batch_clusters),  # Number of clusters
        batch_path  # Current batch path
    ))
    # Launch mgseed jobs
    bjobs = list(map(
        # Mapped function
        lambda cluster_name: run_mgseed(cluster_name, cwd=batch_path, env=env),
        # Input values
        batch_clusters
    ))
    # Check running mgseed.pl jobs
    Bjob.check(bjobs, delay=30)

    # Define set of running jobs
    bjobs = list()
    # Retrieve cluster paths
    cluster_paths = glob.glob(batch_path + '/MGYP*')
    # Debug
    print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
    # Launch pfbuild scripts
    bjobs = list(map(
        # Mapped function
        lambda cluster_path: run_pfbuild(cluster_path, env=env),
        # Input values
        cluster_paths
    ))
    # Check running pfbuild jobs
    Bjob.check(bjobs, delay=30)

    # Run check_uniprot.pl in current batch directory
    out = subprocess.run(
        capture_output=True,
        encoding='utf-8',
        env=env,
        cwd=batch_path,
        args=['check_uniprot.pl']
    )
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
    print('There are {:d} ({:.02f}) remaining clusters: {}'.format(
        # Number of possible MGnifam
        len(kept_clusters),
        # Rate of possible MGnifam clusters
        len(kept_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible MGnifam
        ', '.join(os.path.basename(path) for path in kept_clusters)
    ))
    # Debug
    print('There are {:d} ({:.02f}) clusters in BIAS: {}'.format(
        # Number of possible Pfam
        len(bias_clusters),
        # Rate of possible Pfam clusters
        len(bias_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible Pfam
        ', '.join(os.path.basename(path) for path in bias_clusters)
    ))
    # Debug
    print('There are {:d} ({:.02f}) clusters in Uniprot: {}'.format(
        # Number of possible Pfam
        len(uni_clusters),
        # Rate of possible Pfam clusters
        len(uni_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
        # List all possible Pfam
        ', '.join(os.path.basename(path) for path in uni_clusters)
    ))

    # Delete current batch folder (free space)
    shutil.rmtree(batch_path)

# Define a common multiple sequence alignment transformation pipeline
transform = Compose([
    # Exclude regions outside N- and C- terminal
    OccupancyTrim(
        # Use custom threshold:
        # threshold = lambda x: np.mean(x) - np.std(x),
        threshold=0.5,
        # Set threshold inclusive
        inclusive=True
    ),
    # Exclude sequences with less than half occupancy
    OccupancyFilter(
        # Use custom threshold:
        threshold=0.5,
        # Set threshold inclusive
        inclusive=True
    ),
    # Make non redundant
    MakeNonRedundant(
        # Define redundancy threshold
        threshold=0.8,
        # Set threshold exclusive
        inclusive=False
    )
])
# Initialize summary dictionary
summary = dict()
# Get list of clusters paths
cluster_paths = glob.glob(build_path + '/MGYP*')
# Debug
print('There are {:d} clusters in build: {}'.format(
    len(cluster_paths),
    ', '.join([
        os.path.dirname(os.path.basename(path)) for path in cluster_paths
    ])
))
# For every file in build, run pfbuild again
for cluster_path in cluster_paths:
    # Get cluster name
    cluster_name = os.path.basename(cluster_path)
    cluster_name = os.path.dirname(cluster_name)
    # Initialize summary entry for current cluster
    summary[cluster_name] = {
        # Prettiness score
        'pre_prettiness': None,  # Before trimming
        'post_prettiness': None,  # After trimming
        # Occupancy distribution
        'pre_occupancy': None,  # Before trimming
        'post_occupancy': None,  # After occupancy
        # Conservation distribution
        'pre_conservation': None,  # Before trimming
        'post_conservation': None  # After occupancy
    }
    # Load multiple seed alignment
    seed = MSA().from_aln(cluster_path)
    # Save summary before trimming
    summary[cluster_name] = {**summary[cluster_name], **{
        'pre_prettiness': MSA.prettiness(seed.aln),
        'pre_occupancy': MSA.occupancy(seed.aln)[0],
        'pre_conservation': MSA.conservation(seed.aln)[0]
    }}
    # Execute trim, substitute original seed
    seed = transform(seed)
    # Save summary after trimming
    summary[cluster_name] = {**summary[cluster_name], **{
        'post_prettiness': MSA.prettiness(seed.aln),
        'post_occupancy': MSA.occupancy(seed.aln)[0],
        'post_conservation': MSA.conservation(seed.aln)[0]
    }}
    # Store new SEED alignment
    seed.to_aln(cluster_path)

# Plot summary
fig, axs = plt.subplots(2, 3, figsize=(30, 20), sharex='col', sharey='col')
# Set titles
_ = axs[0, 0].set_title('Pre-trim prettiness')
_ = axs[0, 1].set_title('Pre-trim occupancy')
_ = axs[0, 2].set_title('Pre-trim conservation')
_ = axs[1, 0].set_title('Post-trim prettiness')
_ = axs[1, 1].set_title('Post-trim occupancy')
_ = axs[1, 2].set_title('Post-trim conservation')
# Plot prettiness (pre-trim)
pre_prettiness = [summary[k]['pre_prettiness'] for k in summary.keys()]
_ = axs[0, 0].hist(pre_prettiness, bins=100)
_ = axs[0, 0].axvline(np.mean(pre_prettiness))
# Plot prettiness (post-trim)
post_prettiness = [summary[k]['post_prettiness'] for k in summary.keys()]
_ = axs[1, 0].hist(post_prettiness, bins=100)
_ = axs[1, 0].axvline(np.mean(post_prettiness))
# Plot occupancy distribution (pre-trim)
for k in summary.keys():
    _ = axs[0, 1].hist(
        summary[k]['pre_occupancy'],
        density=True,
        stacked=True,
        bisn=100
    )
# Plot occupancy distribution (post-trim)
for k in summary.keys():
    _ = axs[1, 1].hist(
        summary[k]['post_occupancy'],
        density=True,
        stacked=True,
        bins=100
    )
# Plot conservation distribution (pre-trim)
for k in summary.keys():
    _ = axs[0, 2].hist(
        summary[k]['pre_conservation'],
        density=True,
        stacked=True,
        bins=100
    )
# Plot conservation distribution (post-trim)
for k in summary.keys():
    _ = axs[1, 2].hist(
        summary[k]['post_conservation'],
        density=True,
        stacked=True,
        bins=100
    )
# Save figure to file
_ = plt.savefig(build_path + '/trim.png')
# Close plot
_ = plt.close()
# Delete summary variable (free some space)
del summary

# Debug
print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
# Launch pfbuild jobs
bjobs = list(map(
    # Mapped function
    lambda cluster_path: run_pfbuild(cluster_path, env=env),
    # Input values list
    cluster_paths
))
# Check running pfbuild jobs
Bjob.check(bjobs, delay=30)

# Run check_uniprot.pl in current build directory
out = subprocess.run(
    capture_output=True,
    encoding='utf-8',
    env=env,
    cwd=build_path,
    args=['check_uniprot.pl']
)
# Debug
print('check_uniprot.pl:', out)

# Define kept clusters (MGnifam)
kept_clusters = glob.glob(build_path + '/MGYP*')
# Define discarded clusters
bias_clusters = glob.glob(build_path + '/BIAS/MGYP*')
uni_clusters = glob.glob(build_path + '/Uniprot/MGYP*')
# Debug
print('There are {:d} ({:.02f}) remaining clusters: {}'.format(
    # Number of possible MGnifam
    len(kept_clusters),
    # Rate of possible MGnifam clusters
    len(kept_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
    # List all possible MGnifam
    ', '.join(os.path.basename(path) for path in kept_clusters)
))
# Debug
print('There are {:d} ({:.02f}) clusters in BIAS: {}'.format(
    # Number of possible Pfam
    len(bias_clusters),
    # Rate of possible Pfam clusters
    len(bias_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
    # List all possible Pfam
    ', '.join(os.path.basename(path) for path in bias_clusters)
))
# Debug
print('There are {:d} ({:.02f}) clusters in Uniprot: {}'.format(
    # Number of possible Pfam
    len(uni_clusters),
    # Rate of possible Pfam clusters
    len(uni_clusters) / (len(bias_clusters) + len(uni_clusters) + len(kept_clusters)),
    # List all possible Pfam
    ', '.join(os.path.basename(path) for path in uni_clusters)
))
# Remove clusters in BIAS
shutil.rmtree(build_path + '/BIAS')
# Remove clusters in Uniprot
shutil.rmtree(build_path + '/Uniprot')
