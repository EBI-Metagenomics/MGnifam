# Common dependencies
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import subprocess
import argparse
import shutil
import glob
import time
import json
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


# Run mgseed.pl
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


# Run pfbuild/mfbuild
def run_pfbuild(cluster_path, db='uniprot', withpfmake=True, make_eval=None, env=os.environ.copy()):
    # Make list of arguments
    args = ['pfbuild']
    # Check if 'withpfmake' is enabled
    args += ['-withpfmake'] if withpfmake else []
    # Define database
    args += ['-db', db]
    # Check if 'make eval' is enabled
    args += ['-makeEval', str(make_eval)] if make_eval is not None else []
    # Run pfbuild current cluster directory
    out = subprocess.run(
        capture_output=True,  # Capture console output
        encoding='utf-8',  # Set output encoding
        env=env,  # Set custom environment
        cwd=cluster_path,  # Set directory
        args=args
    )
    # Debug
    print('pfbuild:', out)
    # Get process id as string
    job_id = Bjob.id_from_string(out.stdout)
    # Return new Bjob instance
    return Bjob(id=job_id, status='RUN')


# Wrapper for running mfbuild
def run_mfbuild(cluster_path, env=os.environ.copy()):
    return run_pfbuild(
        cluster_path=cluster_path,
        db='mgnify',
        withpfmake=True,
        make_eval=0.01,
        env=env
    )


# Run check_uniprot.pl
def run_check_uniprot(clusters_dir, env=os.environ.copy()):
    # Run check_uniprot.pl in current batch directory
    out = subprocess.run(
        capture_output=True,
        encoding='utf-8',
        env=env,
        cwd=clusters_dir,
        args=['check_uniprot.pl']
    )
    # Debug
    print('check_uniprot.pl:', out)


# Make batch of clusters
def make_batch_clusters(cluster_names, batch_path, env=os.environ.copy(), log=dict()):
    """
    Firstly, makes a seed alignment for each cluster in the batch by running
    mgseed.pl (discard biased clusters by putting them in BIAS/ direcory).
    Then, makes HMMs by running pfbuild. Finally, discards possible Pfam
    clusters by putting them in Uniprot/ directory.

    Args
    cluster_names (list):   List of cluster names to search in clusters file
    batch_path (str):       Path to current batch directory, where output will
                            be stored
    env (dict)              Environmental variables, required for running
                            mgseed.pl and pfbuild safely
    """
    # Log
    log['beg_time'] = time.time()  # Define iteration start time

    # Make batch directory
    os.mkdir(batch_path)

    # Debug
    print('Running mgseed.pl for all the {} clusters in batch {}'.format(
        len(cluster_names),  # Number of clusters
        batch_path  # Current batch path
    ))
    # Launch mgseed jobs
    bjobs = list(map(
        # Mapped function
        lambda cluster_name: run_mgseed(cluster_name, cwd=batch_path, env=env),
        # Input values
        cluster_names
    ))
    # Check running mgseed.pl jobs
    Bjob.check(bjobs, delay=30)

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
    run_check_uniprot(batch_path, env=env)

    # Log
    log['end_time'] = time.time()  # Set iteation end time
    log['took_time'] = log['end_time'] - log['beg_time']  # Compute time for iteration

    # Define MGnifam (kept) clusters
    mgnifam_clusters = glob.glob(batch_path + '/MGYP*')
    # Define BIAS (discarded) clusters
    bias_clusters = glob.glob(batch_path + '/BIAS/MGYP*')
    # Define Uniprot (discarded) clusters
    uniprot_clusters = glob.glob(batch_path + '/Uniprot/MGYP*')

    # Log
    log['num_mgnifam'] = len(mgnifam_clusters)  # Number of MGnifam clusters
    log['num_bias'] = len(bias_clusters)  # Number of BIAS discarded clusters
    log['num_uniprot'] = len(uniprot_clusters)  # Number of Uniprot discarded clusters

# Initialize log (dictionary)
log = {
    # Begininning time of the script
    'beg_time':  time.time(),
    # End time of the script
    'end_time':  None,
    # How much did it take to execute (seconds)
    'took_time':  0.0,
    # Initialize log for every batch (list), each batch formatted as:
    #   - Batch iteration index
    #   - Start time of iteration
    #   - End time of iteration
    #   - How much time a single batch did take (seconds)
    #   - Number of Mgnifam entries in current batch
    #   - Number of Pfam entries in current batch
    #   - Number of BIAS entries in current batch
    'batch_iter': [],
    # Number of MGnifam entries at the end
    'num_mgnifam': 0,
    # Number of Pfam entries at the end
    'num_uniprot': 0,
    # Number of BIAS entries at the end
    'num_bias': 0,
    # User-defined batch size
    'batch_size': 0,
    # User-defined maximum number of clusters
    'num_clusters': 0
}


# Read arguments
parser = argparse.ArgumentParser(description='Develop a new MGnifam release')
parser.add_argument(
    '-i', '--in_clusters', nargs='+', type=str, required=True,
    help='Path to clusters file, returned by LinClust'
)
parser.add_argument(
    '-o', '--out_dir', type=str, required=True,
    help='Path to directory where outcomes will be stored'
)
parser.add_argument(
    '--num_clusters', type=int, default=None,
    help='Maximum number of clusters to process'
)
parser.add_argument(
    '--batch_size', type=int, default=1000,
    help='Number of clusters which can be processed at the same time'
)
args = parser.parse_args()

# Get output directory path
out_dir = re.sub(r'/$', '', args.out_dir)

# Get input clusters paths
# in_clusters = glob.glob(args.in_clusters)
in_clusters = args.in_clusters
# Get maximum number of cluster to process
num_clusters = args.num_clusters
num_clusters = np.inf if num_clusters is None else int(num_clusters)
# # Get number of clusters
# num_clusters = len(in_clusters)
# Get batch size
batch_size = args.batch_size
# Define build directory
build_path = out_dir + '/build'
# Make build directory
os.mkdir(build_path)

# Update log
log = {**log, **{
    'batch_size': batch_size,
    'num_clusters': num_clusters
}}

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

# Initialize cluster index (points to current cluster name)
cluster_idx = 0
# Initialize list containing cluster names of size batch size
batch_clusters = list()
# Loop through every given input file path
for in_path in in_clusters:
    # Open file reader buffer
    in_file = open(in_path, 'r')
    # Loop through every line of the currently looped input file
    for line in in_file:
        # Check early stopping condition
        # Note: cluster index holds the number of clusters previously looped
        if cluster_idx >= num_clusters:
            break  # Exit inner loop
        # Match cluster name and size in current line
        match = re.search(r'^([a-zA-Z0-9]+)[ \t]+(\d+)', line)
        # Case current file line does not match format, go to next line
        if not match:
            continue  # Skip iteration, go to next line
        # Case current line matches format, store cluster name and size
        cluster_name, cluster_size = str(match.group(1)), int(match.group(2))
        # Add cluster name to current batch
        batch_clusters.append(cluster_name)
        # Update cluster index
        cluster_idx += 1
        # Case the number of cluster names in batch maatches batch size
        if (cluster_idx % batch_size) == 0:
            # Initialize log for current batch
            batch_log = {'batch_idx': cluster_idx // batch_size}
            # Make seed alignment and HMMs for current batch
            _ = make_batch_clusters(
                batch_clusters,  # Batch of cluster names
                # Path to current batch
                batch_path=out_dir+'/batch{:d}'.format(
                    # Compute number of clusters batch
                    cluster_idx // batch_size
                ),
                env=env,  # Environmental variables
                log=batch_log  # Log
            )
            # Save in current log
            log['batch_iter'].append(batch_log)
            # Initialize new empty batch
            batch_clusters = list()
        # Case the number of clusters match the early stopping condition
        if cluster_idx >= num_clusters:
            # Case there is at least one element in new batch
            if (cluster_idx % batch_size) != 0:
                # Initialize log for current batch
                batch_log = {'batch_idx': (cluster_idx // batch_size) + 1}
                # Make seed alignment and HMMs for remaining batch
                _ = make_batch_clusters(
                    batch_clusters,  # Batch of cluster names
                    # Path to current batch
                    batch_path=out_dir+'/batch{:d}'.format(
                        # Compute number of last clusters batch
                        (cluster_idx // batch_size) + 1
                    ),
                    env=env,  # Environmental variables
                    log=batch_log  # Log
                )
                # Save in current log
                log['batch_iter'].append(batch_log)
    # Close current file
    in_file.close()
    # Check early stopping condition, as in inner loop
    if cluster_idx >= num_clusters:
        break  # Exit outer loop

# Loop through folders not discarded from each batch
for cluster_path in glob.glob(out_dir + '/batch*/MGYP*'):
    # Get cluster name
    cluster_name = os.path.basename(cluster_path)
    # cluster_name = os.path.dirname(cluster_name)
    # Move cluster to build folder
    shutil.move(cluster_path, build_path + '/' + cluster_name)
    # Debug
    print('Moved {} to {}'.format(cluster_path, build_path + '/' + cluster_name))

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
        threshold=0.4,
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
    # Number of folders
    len(cluster_paths),
    # List all folder names
    ', '.join([os.path.basename(path) for path in cluster_paths])
))

# Initialize occupancy
pre_occ, post_occ = list(), list()
# Initialize conservation
pre_cnv, post_cnv = list(), list()
# Initialize prettiness
pre_pty, post_pty = list(), list()
# For every file in build, run pfbuild again
for cluster_path in cluster_paths:
    # Get cluster name
    cluster_name = os.path.basename(cluster_path)
    # Load multiple seed alignment from file
    seed = MSA().from_aln(cluster_path + '/SEED')

    # Update occupancy, conervation and prettyness before trimming
    pre_pty.append(MSA.prettiness(seed.aln))
    pre_occ.extend(MSA.occupancy(seed.aln)[0])
    pre_cnv.extend(MSA.conservation(seed.aln)[0])

    # Execute trimming, substitute original seed
    seed = transform(seed)

    # Update occupancy, conervation and prettyness after trimming
    post_pty.append(MSA.prettiness(seed.aln))
    post_occ.extend(MSA.occupancy(seed.aln)[0])
    post_cnv.extend(MSA.conservation(seed.aln)[0])

    # Store new SEED alignment
    seed.to_aln(cluster_path + '/SEED')

# Plot summary
fig, axs = plt.subplots(2, 3, figsize=(30, 15), sharex='col', sharey='col')
# Set titles
_ = axs[0, 0].set_title('Pre-trim prettiness')
_ = axs[0, 1].set_title('Pre-trim occupancy')
_ = axs[0, 2].set_title('Pre-trim conservation')
_ = axs[1, 0].set_title('Post-trim prettiness')
_ = axs[1, 1].set_title('Post-trim occupancy')
_ = axs[1, 2].set_title('Post-trim conservation')
# Plot prettiness (pre-trim)
_ = axs[0, 0].hist(x=pre_pty, bins=100)
_ = axs[0, 0].axvline(np.mean(pre_pty), color='r')
# Plot prettiness (post-trim)
_ = axs[1, 0].hist(x=post_pty, bins=100)
_ = axs[1, 0].axvline(np.mean(post_pty), color='r')
# Plot occupancy distribution (pre-trim)
_ = axs[0, 1].hist(
    x=pre_occ,
    density=True,
    bins=100
)
_ = axs[0, 1].set_xlim(left=0.0, right=1.0)
# Plot occupancy distribution (post-trim)
_ = axs[1, 1].hist(
    x=post_occ,
    density=True,
    bins=100
)
_ = axs[1, 1].set_xlim(left=0.0, right=1.0)
# Plot conservation distribution (pre-trim)
_ = axs[0, 2].hist(
    x=pre_cnv,
    density=True,
    bins=100
)
_ = axs[0, 2].set_xlim(left=0.0)
# Plot conservation distribution (post-trim)
_ = axs[1, 2].hist(
    x=post_cnv,
    density=True,
    bins=100
)
_ = axs[1, 2].set_xlim(left=0.0)
# Save figure to file
_ = plt.savefig(out_dir + '/trim.png')
# Close plot
_ = plt.close()

# Debug
print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
# Launch pfbuild jobs
Bjob.check(delay=30, bjobs=list(map(
    # Mapped function
    lambda cluster_path: run_pfbuild(cluster_path, env=env),
    # Input values list
    cluster_paths
)))

# Run check_uniprot.pl in current build directory
run_check_uniprot(build_path, env=env)

# Retrieve cluster paths
cluster_paths = glob.glob(build_path + '/MGYP*')
# Debug
print('Clusters for mfbuild:\n{}'.format('\n'.join(cluster_paths)))
# Launch mfbuild jobs
Bjob.check(delay=30, bjobs=list(map(
    # Mapped function
    lambda cluster_path: run_mfbuild(cluster_path, env=env),
    # Input values
    cluster_paths
)))

# Log
log['end_time'] = time.time()  # End time of the script
log['took_time'] = log['end_time'] - log['beg_time']  # Time took by the whole script

# Define MGnifam (kept) clusters
mgnifam_clusters = glob.glob(build_path + '/MGYP*')
# Define BIAS (discarded) clusters
bias_clusters = glob.glob(build_path + '/BIAS/MGYP*')
# Define Uniprot (discarded) clusters
uniprot_clusters = glob.glob(build_path + '/Uniprot/MGYP*')

# Log
log['num_mgnifam'] = len(mgnifam_clusters)  # Number of MGnifam clusters
log['num_bias'] = len(bias_clusters)  # Number of BIAS discarded clusters
log['num_uniprot'] = len(uniprot_clusters)  # Number of Uniprot discarded clusters

# Save log to disk
with open(out_dir + '/log.json', 'w') as log_file:
    # Write out json, indented
    json.dump(log, log_file, indent=True)

# # Remove clusters in BIAS
# shutil.rmtree(build_path + '/BIAS')
# # Remove clusters in Uniprot
# shutil.rmtree(build_path + '/Uniprot')
