# Common dependencies
from os import path, environ, mkdir
from glob import glob, iglob
from time import time
import numpy as np
import subprocess
import argparse
import shutil
import glob
import json
import sys
import os
import re

# Updating python path
sys.path.append(path.dirname(path.realpath(__file__) + '/..'))

# Custom dependencies
from src.pipeline.pipeline import Pipeline, Log
from src.utils import benchmark, get_paths
from src.bjob import Bjob
from src.msa import MSA
from src.transform import Transform
from src.transform import Compose
from src.transform import OccupancyTrim
from src.transform import OccupancyFilter


# Define new pipeline
class Release(Pipeline):

    # Constructor
    def __init__(self, delay=30, env=environ.copy()):
        # Set LSF jobs checking delay, in seconds
        self.delay = delay
        # Store environmental variables
        self.env = env

    # Execute pipeline
    def run(self, clusters_paths, out_dir, batch_size=100, max_clusters=None, log_path='', verbose=False):

        # Case output directory does not exist: make one
        if not path.exists(out_dir):
            # Make directpry (could raise filesystem errors)
            mkdir(out_dir)

        # Intiialize timer
        tot_time = time()
        # Initialize logging dictionary
        log = Log(log_path=log_path, log_dict={
            # Input parameters
            'batch_size': batch_size,
            'max_clusters': max_clusters,
            'out_dir': out_dir,
            # pfbuild computational time
            'pfbuild_time': 0.0,
            # mfbuild computational time
            'mfbuild_time': 0.0,
            # Total execution time
            'tot_time': 0.0
        })

        # Define batches iterator
        batch_iter = self.iter_clusters(clusters_paths, batch_size=batch_size, max_clusters=max_clusters)
        # Loop through each batch
        for batch_index, cluster_names in batch_iter:
            # Define batch directory path
            batch_path = path.join(out_dir, 'batch_{:d}'.format(batch_index))
            # Run batch clusters maker
            self.run_batch(
                cluster_names=cluster_names,  # List of cluster names
                batch_path=batch_path  # Clusters output directory
            )

        # Define build directory
        build_path = path.join(out_dir, 'build')
        # Make build directory
        mkdir(build_path)

        # Reset paths to clusters
        clusters_paths = list()
        # Copy all clusters kept in batch directories to build directories
        for cluster_batch_path in iglob(path.join(out_dir, 'batch*', 'MGYP*')):
            # Get cluster name
            cluster_name = path.basename(cluster_batch_path)
            # Define new cluster path
            cluster_build_path = path.join(build_path, cluster_name)
            # Copy from batch path to build path
            shutil.copytree(cluster_batch_path, cluster_build_path)
            # Save cluster build path to
            clusters_paths.append(cluster_build_path)

        # Define SEED multiple sequence alignments transformations
        transform = Compose([
            # Exclude regions outside N- and C- terminal
            OccupancyTrim(threshold=0.4, inclusive=True),
            # Exclude sequences with less than half occupancy
            OccupancyFilter(threshold=0.5, inclusive=True)
        ])

        # Intialize path for discarded clusters
        noaln_path = os.path.join(build_path, 'NOALN')
        # Make directory
        mkdir(noaln_path)
        # Loop through each SEED alignment
        for cluster_path in clusters_paths:
            # Trim cluster's SEED alignment
            self.trim_seed(
                cluster_path=cluster_path,
                noaln_path=noaln_path,
                transform=transform
            )

        # Run HMMs against UniProt
        _, took_time = benchmark(fn=Bjob.check, delay=self.delay, bjobs=list(map(
            # Mapped function
            lambda cluster_path: self.pfbuild(cluster_path),
            # Input values
            clusters_paths
        )))
        # Store execution time
        log({'pfbuild_time': round(took_time, 2)})

        # Check uniprot overlappings in build path
        self.check_uniprot(build_path)

        # Reset paths to clusters
        clusters_paths = glob(path.join(build_path, 'MGYP*'))
        # Run HMMs against MGnify
        _, took_time = benchmark(fn=Bjob.check, delay=self.delay, bjobs=list(map(
            # Mapped function
            lambda cluster_path: self.mfbuild(cluster_path),
            # Input values
            clusters_paths
        )))
        # Store execution time
        log({'mfbuild_time': round(took_time, 2)})

        # Update log
        log({'tot_time': round(time() - tot_time, 2)})

    @staticmethod
    def iter_clusters(clusters_paths, batch_size=1000, max_clusters=None):
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
        # Loop through each linclust input path
        for cluster_path in get_paths(clusters_paths):
            # Open current linclust path
            with open(cluster_path, 'r') as cluster_file:
                # Loop through each line in current linclust file path
                for cluster_line in cluster_file:
                    # Check if current line is matches expected format
                    match = re.search(r'^(\S+)\s+(\d+)', cluster_line)
                    # Case line does not match expected format: skip
                    if not match: continue
                    # Otherwise, retrieve cluster name and cluster size
                    cluster_name = str(match.group(1))
                    cluster_size = str(match.group(2))
                    # Add current cluster name to batch
                    batch_clusters.append(cluster_name)
                    # Update current cluster index
                    cluster_index += 1
                    # Define batch index
                    batch_index = (cluster_index - 1) // batch_size
                    # Case cluster index has reached maximum size, exit loop
                    if (max_clusters is not None) and (cluster_index >= max_clusters):
                        break
                    # Check if cluster index has not reached batch size
                    elif (cluster_index % batch_size) != 0:
                        # Go to next iteration
                        continue
                    # Case cluster size has reached batch size
                    else:
                        # Yield batch index and list of clusters
                        yield batch_index, batch_clusters
                        # Reset batch of clusters
                        batch_clusters = list()
                # Case cluster index has reached maximum size, exit loop
                if (max_clusters is not None) and (cluster_index >= max_clusters):
                    break
        # In case we reached this point, check for non returned clusters
        if len(batch_clusters) > 0:
            # Yield last batch
            yield batch_index, batch_clusters

    # Run batch of clusters
    def run_batch(self, cluster_names, batch_path, verbose=False):
        """Make batch of HMMs

        Firstly, makes a seed alignment for each cluster in the batch by running
        mgseed.pl (discard biased clusters by putting them in BIAS/ direcory).

        Afterwards, makes HMMs by running pfbuild and discards possible Pfam
        clusters by putting them in Uniprot/ directory.

        Args
        cluster_names (list):   List of cluster names to search in clusters file
        batch_path (str):       Path to current batch directory, where output
                                will be stored
        env (dict)              Environmental variables, required for running
                                mgseed.pl and pfbuild safely

        Return
        (dict)                  Log dictionary containing execution times

        Raise
        (OSError)               In case of permission or filesystem issues
        """

        # # Debug
        # print('Running mgseed.pl for {:d} clusters in batch {:s}'.format(
        #     len(cluster_names),  # Number of clusters
        #     batch_path  # Current batch path
        # ))

        # Make batch directory (could raise error)
        mkdir(batch_path)
        # Start timer
        tot_time = time()
        # Define log path
        log_path = os.path.join(batch_path, 'batch.json')
        # Initialize log
        log = Log(log_path=log_path, log_dict={
            'pfbuild_time': 0.0,
            'mfbuild_time': 0.0,
            'tot_time': 0.0
        })

        # Launch and check mgseed jobs
        _, took_time = benchmark(fn=Bjob.check, delay=self.delay, bjobs=list(map(
            # Mapped function
            lambda cluster_name: self.mgseed(
                cluster_name,  # Current cluster name
                cwd=batch_path  # Current batch path
            ),
            # Input values
            cluster_names
        )))
        # Update log
        log({'mgseed_time': round(took_time, 2)})

        # Retrieve cluster paths
        cluster_paths = glob.glob(batch_path + '/MGYP*')
        # Launch pfbuild scripts
        _, took_time = benchmark(fn=Bjob.check, delay=self.delay, bjobs=list(map(
            # Mapped function
            lambda cluster_path: self.pfbuild(
                cluster_path
            ),
            # Input values
            cluster_paths
        )))
        # Update log
        log({'pfbuild_time': round(took_time, 2)})

        # Check uniprot overlappings
        self.check_uniprot(batch_path)

        # Update log
        log({'tot_time': round(tot_time - time(), 2)})

    # Run mgseed.pl
    def mgseed(self, cluster_name, cwd='./'):
        # Run mgseed.pl in current batch directory
        ran = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            env=self.env,
            cwd=cwd,  # Set directory
            args=['mgseed.pl', '-cluster', cluster_name]
        )
        # Get process id as string
        job_id = Bjob.id_from_string(ran.stdout)
        # Return new Bjob instance
        return Bjob(id=job_id, status='RUN')

    # Run pfbuild/mfbuild
    def pfbuild(self, cluster_path, db='uniprot', withpfmake=True, make_eval=None):
        # Make list of arguments
        args = ['pfbuild']
        # Check if 'withpfmake' is enabled
        args += ['-withpfmake'] if withpfmake else []
        # Define database
        args += ['-db', db]
        # Check if 'make eval' is enabled
        args += ['-makeEval', str(make_eval)] if make_eval is not None else []
        # Run pfbuild current cluster directory
        ran = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set custom environment
            cwd=cluster_path,  # Set directory
            args=args
        )
        # Get process id as string
        job_id = Bjob.id_from_string(ran.stdout)
        # Return new Bjob instance
        return Bjob(id=job_id, status='RUN')

    # Run pfbuild with mfbuild parameters as default
    def mfbuild(self, cluster_path, db='mgnify', withpfmake=True, make_eval=0.01):
        # Just pass parameters to pfbuild script
        return self.pfbuild(*args, **kwargs)

    # Run check_uniprot.pl
    def check_uniprot(self, clusters_dir):
        # Run check_uniprot.pl in current batch directory
        ran = subprocess.run(
            capture_output=True,
            encoding='utf-8',
            env=self.env,
            cwd=clusters_dir,
            args=['check_uniprot.pl']
        )

    # Trim SEED msa
    def trim_seed(self, cluster_path, noaln_path, transform):
        """Trim alignments and plot results
        For each cluster in <clusters_paths>, take SEED alignment and rename it as
        SEED_original. Afterwards, create a new SEED alignment by automatically
        trimming the original one. Automatic trimming is done by feeding alignments
        to transform pipeline.

        Args
        cluster_path (list(str))        Path to cluster which must be trimmed
        noaln_path (list(str))          Path where empty trimmed alignments
                                        must be moved (discarded)
        transform (transform.Transform) Transformer which will be fed with original
                                        SEED alignment and produce the trimmed one
        out_path (str)                  Path where to store generic outputs such as
                                        trimming summary plot
        """
        # Get cluster name from cluster path
        cluster_name = path.basename(cluster_path)
        # Rename original cluster SEED alignment
        shutil.move(
            # Move from
            path.join(cluster_path, 'SEED'),
            # Move to
            path.join(cluster_path, 'SEED_raw')
        )
        # Load raw SEED multiple sequence alignment from file
        seed = MSA.from_aln(path.join(cluster_path, 'SEED_raw'))
        # Go through multiple sequence alignment transformation pipeline
        seed = transform(seed)
        # Case trimmed alignment is empty
        if seed.is_empty():
            # Discard current cluster
            shutil.move(cluster_path, os.path.join(noaln_path, cluster_name))
        # Otherwise
        else:
            # Store new SEED multiple sequence alignment to file
            seed.to_aln(path.join(cluster_path, 'SEED'))


# # Run mgseed.pl
# def run_mgseed(cluster_name, cwd='./', env=os.environ.copy()):
#     # Run mgseed.pl in current batch directory
#     out = subprocess.run(
#         capture_output=True,  # Capture console output
#         encoding='utf-8',  # Set output encoding
#         env=env,
#         cwd=cwd,  # Set directory
#         args=['mgseed.pl', '-cluster', cluster_name]
#     )
#     # Debug
#     print('mgseed.pl:', out)
#     # Get process id as string
#     job_id = Bjob.id_from_string(out.stdout)
#     # Return new Bjob instance
#     return Bjob(id=job_id, status='RUN')


# # Run pfbuild/mfbuild
# def run_pfbuild(cluster_path, db='uniprot', withpfmake=True, make_eval=None, env=os.environ.copy()):
#     # Make list of arguments
#     args = ['pfbuild']
#     # Check if 'withpfmake' is enabled
#     args += ['-withpfmake'] if withpfmake else []
#     # Define database
#     args += ['-db', db]
#     # Check if 'make eval' is enabled
#     args += ['-makeEval', str(make_eval)] if make_eval is not None else []
#     # Run pfbuild current cluster directory
#     out = subprocess.run(
#         capture_output=True,  # Capture console output
#         encoding='utf-8',  # Set output encoding
#         env=env,  # Set custom environment
#         cwd=cluster_path,  # Set directory
#         args=args
#     )
#     # Debug
#     print('pfbuild:', out)
#     # Get process id as string
#     job_id = Bjob.id_from_string(out.stdout)
#     # Return new Bjob instance
#     return Bjob(id=job_id, status='RUN')


# # Wrapper for running mfbuild
# def run_mfbuild(cluster_path, env=os.environ.copy()):
#     return run_pfbuild(
#         cluster_path=cluster_path,
#         db='mgnify',
#         withpfmake=True,
#         make_eval=0.01,
#         env=env
#     )


# # Run check_uniprot.pl
# def run_check_uniprot(clusters_dir, env=os.environ.copy()):
#     # Run check_uniprot.pl in current batch directory
#     out = subprocess.run(
#         capture_output=True,
#         encoding='utf-8',
#         env=env,
#         cwd=clusters_dir,
#         args=['check_uniprot.pl']
#     )
#     # Debug
#     print('check_uniprot.pl:', out)


# # Make batch of clusters
# def make_batch_clusters(cluster_names, batch_path, env=os.environ.copy(), log=dict()):
#     """
#     Firstly, makes a seed alignment for each cluster in the batch by running
#     mgseed.pl (discard biased clusters by putting them in BIAS/ direcory).
#     Then, makes HMMs by running pfbuild. Finally, discards possible Pfam
#     clusters by putting them in Uniprot/ directory.
#
#     Args
#     cluster_names (list):   List of cluster names to search in clusters file
#     batch_path (str):       Path to current batch directory, where output will
#                             be stored
#     env (dict)              Environmental variables, required for running
#                             mgseed.pl and pfbuild safely
#     """
#     # Log
#     log['beg_time'] = time.time()  # Define iteration start time
#
#     # Make batch directory
#     mkdir(batch_path)
#
#     # Debug
#     print('Running mgseed.pl for all the {} clusters in batch {}'.format(
#         len(cluster_names),  # Number of clusters
#         batch_path  # Current batch path
#     ))
#     # Launch mgseed jobs
#     bjobs = list(map(
#         # Mapped function
#         lambda cluster_name: run_mgseed(cluster_name, cwd=batch_path, env=env),
#         # Input values
#         cluster_names
#     ))
#     # Check running mgseed.pl jobs
#     Bjob.check(bjobs, delay=30)
#
#     # Retrieve cluster paths
#     cluster_paths = glob.glob(batch_path + '/MGYP*')
#     # Debug
#     print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
#     # Launch pfbuild scripts
#     bjobs = list(map(
#         # Mapped function
#         lambda cluster_path: run_pfbuild(cluster_path, env=env),
#         # Input values
#         cluster_paths
#     ))
#     # Check running pfbuild jobs
#     Bjob.check(bjobs, delay=30)
#
#     # Run check_uniprot.pl in current batch directory
#     run_check_uniprot(batch_path, env=env)
#
#     # Log
#     log['end_time'] = time.time()  # Set iteation end time
#     log['took_time'] = log['end_time'] - log['beg_time']  # Compute time for iteration
#
#     # Define MGnifam (kept) clusters
#     mgnifam_clusters = glob.glob(batch_path + '/MGYP*')
#     # Define BIAS (discarded) clusters
#     bias_clusters = glob.glob(batch_path + '/BIAS/MGYP*')
#     # Define Uniprot (discarded) clusters
#     uniprot_clusters = glob.glob(batch_path + '/Uniprot/MGYP*')
#
#     # Log
#     log['num_mgnifam'] = len(mgnifam_clusters)  # Number of MGnifam clusters
#     log['num_bias'] = len(bias_clusters)  # Number of BIAS discarded clusters
#     log['num_uniprot'] = len(uniprot_clusters)  # Number of Uniprot discarded clusters


# # Trim SEED alignments
# def trim_alignments(clusters_paths, transform, out_path):
#     """Trim alignments and plot results
#     For each cluster in <clusters_paths>, take SEED alignment and rename it as
#     SEED_original. Afterwards, create a new SEED alignment by automatically
#     trimming the original one. Automatic trimming is done by feeding alignments
#     to transform pipeline.
#
#     Args
#     clusters_paths (list(str))      List of clusters directories, named after
#                                     cluster name itself
#     transform (transform.Transform) Transformer which will be fed with original
#                                     SEED alignment and produce the trimmed one
#     out_path (str)                  Path where to store generic outputs such as
#                                     trimming summary plot
#     """
#     # Initialize occupancy
#     pre_occ, post_occ = list(), list()
#     # Initialize conservation
#     pre_cnv, post_cnv = list(), list()
#     # Initialize prettiness
#     pre_pty, post_pty = list(), list()
#     # For every file in build, run pfbuild again
#     for cluster_path in cluster_paths:
#         # Get cluster name
#         cluster_name = os.path.basename(cluster_path)
#         # Rename original cluster SEED alignment
#         shutil.move(cluster_path + '/SEED', cluster_path + '/SEED_raw')
#         # Load multiple seed alignment from file
#         seed = MSA.from_aln(cluster_path + '/SEED_original')
#
#         # Initialize pre-trimming plot
#         fig = plt.figure(constrained_layout=True, figsize=(20, 10))
#         # Define axes grid
#         grid = fig.add_gridspec(2, 4)
#         # Add main title
#         fig.suptitle('Pre-trimming multiple sequence alignment (MSA)')
#         # Initialize axis
#         axs = [
#             fig.add_subplot(grid[0, 0]),
#             fig.add_subplot(grid[0, 1]),
#             fig.add_subplot(grid[0, 2:]),
#             fig.add_subplot(grid[1, 0]),
#             fig.add_subplot(grid[1, 1]),
#             fig.add_subplot(grid[1, 2:])
#         ]
#         # Make plots in axes
#         Transform.plot_results(msa=seed, axs=axs)
#         # Make plot
#         plt.savefig(cluster_path + '/pre_trim.png')
#         plt.close()
#
#         # Update occupancy, conervation and prettyness before trimming
#         pre_pty.append(MSA.prettiness(seed.aln))
#         pre_occ.extend(MSA.occupancy(seed.aln)[0])
#         pre_cnv.extend(MSA.conservation(seed.aln)[0])
#
#         # Execute trimming, substitute original seed
#         seed = transform(seed)
#
#         # Initialize pre-trimming plot
#         fig = plt.figure(constrained_layout=True, figsize=(20, 10))
#         # Define axes grid
#         grid = fig.add_gridspec(2, 4)
#         # Add main title
#         fig.suptitle('Post-trimming multiple sequence alignment (MSA)')
#         # Initialize axis
#         axs = [
#             fig.add_subplot(grid[0, 0]),
#             fig.add_subplot(grid[0, 1]),
#             fig.add_subplot(grid[0, 2:]),
#             fig.add_subplot(grid[1, 0]),
#             fig.add_subplot(grid[1, 1]),
#             fig.add_subplot(grid[1, 2:])
#         ]
#         # Make plots in axes
#         Transform.plot_results(msa=seed, axs=axs)
#         # Make plot
#         plt.savefig(cluster_path + '/post_trim.png')
#         plt.close()
#
#         # Update occupancy, conervation and prettyness after trimming
#         post_pty.append(MSA.prettiness(seed.aln))
#         post_occ.extend(MSA.occupancy(seed.aln)[0])
#         post_cnv.extend(MSA.conservation(seed.aln)[0])
#
#         # Store new SEED alignment
#         seed.to_aln(cluster_path + '/SEED')
#
#     # Plot summary
#     fig, axs = plt.subplots(2, 3, figsize=(30, 15), sharex='col', sharey='col')
#     # Set titles
#     _ = axs[0, 0].set_title('Pre-trim prettiness')
#     _ = axs[0, 1].set_title('Pre-trim occupancy')
#     _ = axs[0, 2].set_title('Pre-trim conservation')
#     _ = axs[1, 0].set_title('Post-trim prettiness')
#     _ = axs[1, 1].set_title('Post-trim occupancy')
#     _ = axs[1, 2].set_title('Post-trim conservation')
#     # Plot prettiness (pre-trim)
#     _ = axs[0, 0].hist(x=pre_pty, bins=100)
#     _ = axs[0, 0].axvline(np.mean(pre_pty), color='r')
#     # Plot prettiness (post-trim)
#     _ = axs[1, 0].hist(x=post_pty, bins=100)
#     _ = axs[1, 0].axvline(np.mean(post_pty), color='r')
#     # Plot occupancy distribution (pre-trim)
#     _ = axs[0, 1].hist(
#         x=pre_occ,
#         density=True,
#         bins=100
#     )
#     _ = axs[0, 1].set_xlim(left=0.0, right=1.0)
#     # Plot occupancy distribution (post-trim)
#     _ = axs[1, 1].hist(
#         x=post_occ,
#         density=True,
#         bins=100
#     )
#     _ = axs[1, 1].set_xlim(left=0.0, right=1.0)
#     # Plot conservation distribution (pre-trim)
#     _ = axs[0, 2].hist(
#         x=pre_cnv,
#         density=True,
#         bins=100
#     )
#     _ = axs[0, 2].set_xlim(left=0.0)
#     # Plot conservation distribution (post-trim)
#     _ = axs[1, 2].hist(
#         x=post_cnv,
#         density=True,
#         bins=100
#     )
#     _ = axs[1, 2].set_xlim(left=0.0)
#     # Save figure to file
#     _ = plt.savefig(out_path + '/trim.png')
#     # Close plot
#     _ = plt.close()


# Main
if __name__ == '__main__':

    # Read arguments
    parser = argparse.ArgumentParser(description='Develop a new MGnifam release')
    parser.add_argument(
        '-i', '--in_clusters', nargs='+', type=str, required=True,
        help='Path to clusters file, returned by LinClust'
    )
    parser.add_argument(
        '-o', '--out_dir', type=str, default='./',
        help='Path to directory where outcomes will be stored'
    )
    parser.add_argument(
        '--max_clusters', type=int, required=False,
        help='Maximum number of clusters to process'
    )
    parser.add_argument(
        '--batch_size', type=int, default=100,
        help='Number of clusters which can be processed at the same time'
    )
    parser.add_argument(
        '--log_path', type=str, default='',
        help='Path where to store log'
    )
    args = parser.parse_args()

    # Get output directory path
    out_dir = re.sub(r'/$', '', args.out_dir)
    # Get input clusters paths
    in_clusters = args.in_clusters
    # Get maximum number of cluster to process
    max_clusters = args.max_clusters
    # Get batch size
    batch_size = args.batch_size
    # Get path to log file
    log_path = args.log_path

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

    # Instantiate new release pipeline
    pipeline = Release(delay=30, env=env)

    # Run new pipeline
    pipeline(
        # Input .tsv file path
        clusters_paths=in_clusters,
        # Output directory
        out_dir=out_dir,
        # Set batch size,
        batch_size=batch_size,
        # Set maximum number of clusters
        max_clusters=max_clusters,
        # Path to log file
        log_path=log_path,
        # Verbose output (not working)
        verbose=True
    )

    # # Initialize cluster index (points to current cluster name)
    # cluster_idx = 0
    # # Initialize list containing cluster names of size batch size
    # batch_clusters = list()
    # # Loop through every given input file path
    # for in_path in in_clusters:
    #     # Open file reader buffer
    #     in_file = open(in_path, 'r')
    #     # Loop through every line of the currently looped input file
    #     for line in in_file:
    #         # Check early stopping condition
    #         # Note: cluster index holds the number of clusters previously looped
    #         if cluster_idx >= num_clusters:
    #             break  # Exit inner loop
    #         # Match cluster name and size in current line
    #         match = re.search(r'^([a-zA-Z0-9]+)[ \t]+(\d+)', line)
    #         # Case current file line does not match format, go to next line
    #         if not match:
    #             continue  # Skip iteration, go to next line
    #         # Case current line matches format, store cluster name and size
    #         cluster_name, cluster_size = str(match.group(1)), int(match.group(2))
    #         # Add cluster name to current batch
    #         batch_clusters.append(cluster_name)
    #         # Update cluster index
    #         cluster_idx += 1
    #         # Case the number of cluster names in batch maatches batch size
    #         if (cluster_idx % batch_size) == 0:
    #             # Initialize log for current batch
    #             batch_log = {'batch_idx': cluster_idx // batch_size}
    #             # Make seed alignment and HMMs for current batch
    #             _ = make_batch_clusters(
    #                 batch_clusters,  # Batch of cluster names
    #                 # Path to current batch
    #                 batch_path=out_dir+'/batch{:d}'.format(
    #                     # Compute number of clusters batch
    #                     cluster_idx // batch_size
    #                 ),
    #                 env=env,  # Environmental variables
    #                 log=batch_log  # Log
    #             )
    #             # Save in current log
    #             log['batch_iter'].append(batch_log)
    #             # Initialize new empty batch
    #             batch_clusters = list()
    #         # Case the number of clusters match the early stopping condition
    #         if cluster_idx >= num_clusters:
    #             # Case there is at least one element in new batch
    #             if (cluster_idx % batch_size) != 0:
    #                 # Initialize log for current batch
    #                 batch_log = {'batch_idx': (cluster_idx // batch_size) + 1}
    #                 # Make seed alignment and HMMs for remaining batch
    #                 _ = make_batch_clusters(
    #                     batch_clusters,  # Batch of cluster names
    #                     # Path to current batch
    #                     batch_path=out_dir+'/batch{:d}'.format(
    #                         # Compute number of last clusters batch
    #                         (cluster_idx // batch_size) + 1
    #                     ),
    #                     env=env,  # Environmental variables
    #                     log=batch_log  # Log
    #                 )
    #                 # Save in current log
    #                 log['batch_iter'].append(batch_log)
    #     # Close current file
    #     in_file.close()
    #     # Check early stopping condition, as in inner loop
    #     if cluster_idx >= num_clusters:
    #         break  # Exit outer loop
    #
    # # Loop through folders not discarded from each batch
    # for cluster_path in glob.glob(out_dir + '/batch*/MGYP*'):
    #     # Get cluster name
    #     cluster_name = os.path.basename(cluster_path)
    #     # Move cluster to build folder
    #     shutil.copy(cluster_path, build_path + '/' + cluster_name)
    #     # Debug
    #     print('Moved {} to {}'.format(cluster_path, build_path + '/' + cluster_name))
    #
    # # Define a common multiple sequence alignment transformation pipeline
    # transform = Compose([
    #     # Exclude regions outside N- and C- terminal
    #     OccupancyTrim(threshold=0.4, inclusive=True),
    #     # Exclude sequences with less than half occupancy
    #     OccupancyFilter(threshold=0.5, inclusive=True)
    # ])
    # # Get list of clusters paths
    # cluster_paths = glob.glob(build_path + '/MGYP*')
    # # Debug
    # print('There are {:d} clusters in build: {}'.format(
    #     # Number of folders
    #     len(cluster_paths),
    #     # List all folder names
    #     ', '.join([os.path.basename(path) for path in cluster_paths])
    # ))
    # # Make trim
    # trim_alignments(cluster_paths, transform, out_dir)
    #
    # # Debug
    # print('Clusters for pfbuild:\n{}'.format('\n'.join(cluster_paths)))
    # # Initialize timers
    #
    # # Launch pfbuild jobs
    # Bjob.check(delay=30, bjobs=list(map(
    #     # Mapped function
    #     lambda cluster_path: run_pfbuild(cluster_path, env=env),
    #     # Input values list
    #     cluster_paths
    # )))
    #
    # # Run check_uniprot.pl in current build directory
    # run_check_uniprot(build_path, env=env)
    #
    # # Retrieve cluster paths
    # cluster_paths = glob.glob(build_path + '/MGYP*')
    # # Debug
    # print('Clusters for mfbuild:\n{}'.format('\n'.join(cluster_paths)))
    # # Launch mfbuild jobs
    # Bjob.check(delay=30, bjobs=list(map(
    #     # Mapped function
    #     lambda cluster_path: run_mfbuild(cluster_path, env=env),
    #     # Input values
    #     cluster_paths
    # )))
    #
    # # Log
    # log['end_time'] = time.time()  # End time of the script
    # log['took_time'] = log['end_time'] - log['beg_time']  # Time took by the whole script
    #
    # # Define MGnifam (kept) clusters
    # mgnifam_clusters = glob.glob(build_path + '/MGYP*')
    # # Define BIAS (discarded) clusters
    # bias_clusters = glob.glob(build_path + '/BIAS/MGYP*')
    # # Define Uniprot (discarded) clusters
    # uniprot_clusters = glob.glob(build_path + '/Uniprot/MGYP*')
    #
    # # Log
    # log['num_mgnifam'] = len(mgnifam_clusters)  # Number of MGnifam clusters
    # log['num_bias'] = len(bias_clusters)  # Number of BIAS discarded clusters
    # log['num_uniprot'] = len(uniprot_clusters)  # Number of Uniprot discarded clusters
    #
    # # Save log to disk
    # with open(out_dir + '/log.json', 'w') as log_file:
    #     # Write out json, indented
    #     json.dump(log, log_file, indent=True)
    #
    # # # Remove clusters in BIAS
    # # shutil.rmtree(build_path + '/BIAS')
    # # # Remove clusters in Uniprot
    # # shutil.rmtree(build_path + '/Uniprot')
