# Common dependencies
from glob import glob, iglob
from time import time
import numpy as np
import subprocess
import argparse
import shutil
import json
import sys
import os
import re

# Updating python path
sys.path.append(os.path.dirname(os.path.realpath(__file__) + '/..'))

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
class OldRelease(Pipeline):

    # Constructor
    def __init__(self, delay=30, env=os.environ.copy()):
        # Set LSF jobs checking delay, in seconds
        self.delay = delay
        # Store environmental variables
        self.env = env
        # Define trimming
        self.seed_transform = Compose([
            # Exclude regions outside N- and C- terminal
            OccupancyTrim(threshold=0.4, inclusive=True),
            # Exclude sequences with less than half occupancy
            OccupancyFilter(threshold=0.5, inclusive=True)
        ])

    # Execute pipeline
    def run(self, in_path, out_path, author_name, batch_size=100, max_clusters=None, verbose=False):
        """Run MGnifam pipeline

        Takes as input a file containing a spaces separated table, whose first
        column contains clusters names.

        Then, makes a first run for each batch of clusters of given batch size,
        until given maximum number of clusters is reached, if any.

        From clusters remaining after the first iteration, checks the remaining
        ones, discards some and makes MGnifam entries from the others.

        Args
        in_path (str)       Path to input file[s]
        out_path (str)      Path to output directory
        author_name (str)   Name of the authore making the new release
        batch_size (int)    How many clusters to run at each iteration
        max_clusters (int)  Maximum number of clusters to go through
        verbose (bool)      Wether to show verbose output or not

        Return
        (float)             Execution time of the pipeline
        (set)               Set of valid cluster names
        """
        # Intialize timers
        time_beg, time_end = time(), 0.0

        # Verbose log
        if verbose:
            print('Making new MGnifam release', end=' ')
            print('at {:s}'.format(out_path))

        # Case output directory does not exist: make one
        if not os.path.exists(out_path):
            # Make output directory (could raise filesystem errors)
            os.mkdir(out_path)

        # Define build directory
        build_path = os.path.join(out_path, 'build')
        # CHeck if directory already exists
        if os.path.exists(build_path):
            # Make build directory
            os.mkdir(build_path)

        # Define batch directory
        batch_path = os.path.join(out_path, 'batch')
        # Case batch directory des not exist
        if not os.path.exists(batch_path):
            # Make one
            os.mkdir(batch_path)

        # Initialize log path
        log_path = os.path.join(out_path, 'release.json')

        # Initialize log
        log = Log(log_path=log_path, log_dict={
            'out_path': out_path,
            'author_name': author_name,
            'batch_size': batch_size,
            'max_clusters': max_clusters,
            'batch_path': batch_path,
            'build_path': build_path
        })

        # Define batches iterator
        batch_iter = self.iter_clusters(in_path, batch_size=batch_size, max_clusters=max_clusters)
        # Initialize batch index
        batch_index = 0
        # Loop through each batch
        for batch_index, cluster_names in batch_iter:
            # Define current batch directory path
            curr_path = os.path.join(batch_path, str(batch_index))
            # Make set out of cluster names list
            cluster_names = set(cluster_names)
            # Run batch
            self.make_batch(
                cluster_names=cluster_names,
                batch_path=curr_path,
                build_path=build_path,
                verbose=verbose
            )

        # Update timer
        batch_end = time()
        batch_tot = batch_end - time_beg

        # Verbose log
        if verbose:
            print('Batches completed: {:d} iterations'.format(batch_index + 1), end=' ')
            print('done in {:.0f} seconds'.format(batch_tot))

        # Update log
        log({
            'num_batches': batch_index + 1,
            'batches_time': batch_tot
        })

        # Make build
        time_run, cluster_names = self.make_build(
            build_path=build_path,
            author_name=author_name,
            verbose=verbose
        )

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('Build done in {:.0f} seconds'.format(time_run))
            print('Release done in {:.0f} seconds:'.format(time_tot), end=' ')
            print('there are {:d} resulting clusters'.format(len(cluster_names)))

        # Update log
        log({
            'build_time': time_run,
            'release_time': time_tot,
            'release_size': len(cluster_names)
        })

        # Return time taken and resulting cluster names
        return time_tot, cluster_names

    @staticmethod
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

    # Run a single batch
    def make_batch(self, cluster_names, batch_path, build_path, verbose=False):
        """Make batch of HMMs

        First, generates SEED alignment taking sequence out of MGnifam dataset,
        checking that their compositional bias threshold is under a certain
        threshold, otherwise discarding them into BIAS/ folder.

        Then, build HMMs out from the remaining clusters. Built HMMs are checked
        against UniProt and if at least one match is found, then cluster is
        discarded by moving it into UniProt/ folder.

        Remaining cluster folders are being copied to build/ folder, to undergo
        another pipeline which will lead to new MGnifam entries.

        Args
        cluster_names (list)        List of cluster names, which will be taken
                                    out from MGnifam dataset
        batch_path (str)            Path to current batch of clusters
        build_path (str)            Path where to move succesfull clusters
        verbose (bool)              Wether to show verbose log or not
        """

        # Start timer
        time_beg, time_end = time(), 0.0

        # Get batch name
        batch_name = os.path.basename(batch_path)

        # Ensure cluster names is a set
        cluster_names = set(cluster_names)

        # Verbose log
        if verbose:
            print('Making batch {:s}'.format(batch_name), end=' ')
            print('of size {:d}'.format(len(cluster_names)))

        # Case batch directory does not exist
        if not os.path.exists(batch_path):
            # Make new directory
            os.mkdir(batch_path)

        # # Define NOALN path
        # noaln_path = os.path.join(batch_path, 'NOALN')
        # # Check if NOALN folder exists
        # if not os.path.exists(noaln_path):
        #     # Try to make new NOALN folder
        #     os.mkdir(noaln_path)

        # Define DISCARD path
        discard_path = os.path.join(batch_path, 'DISCARD')
        # Check if DISCARDED path exist
        if not os.path.exists(discard_path):
            # Try to make new DISCARD folder
            os.mkdir(discard_path)

        # Define log path
        log_path = os.path.join(batch_path, 'batch.json')
        # Initialize log
        log = Log(log_path=log_path, log_dict={
            'batch_name': batch_name,
            'batch_path': batch_path,
            'batch_size': len(cluster_names)
        })

        # Make SEED alignments
        time_run, cluster_names = self.make_seed_alignments(
            clusters_path=batch_path,
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Update log
        log({
            'mgseed_time': time_run,
            'mgseed_size': len(cluster_names)
        })

        # Make HMMs and search them against UniProt
        time_run, cluster_names = self.make_hmm_models(
            clusters_path=batch_path,
            cluster_names=cluster_names,
            db='uniprot',
            withpfmake=True,
            make_eval=None,
            verbose=verbose
        )

        # Update log
        log({
            'pfbuild_time': time_run,
            'pfbuild_size': len(cluster_names)
        })

        # Check HMMs for UniProt insersection
        time_run, cluster_names = self.check_hmm_uniprot(
            clusters_path=batch_path,
            discard_path=discard_path,
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Update log
        log({
            'uniprot_check_time': time_run,
            'uniprot_check_size': len(cluster_names)
        })

        # Define new clusters iterator
        clusters_iter = iglob(os.path.join(batch_path, 'MGYP*'))
        # Copy remaining clusters to build path
        for cluster_path in clusters_iter:
            # Get cluster name
            cluster_name = os.path.basename(cluster_path)
            # Set target path to build directory
            target_path = os.path.join(build_path, cluster_name)
            # Case cluster is not in direcory
            if not os.path.exists(cluster_path):
                continue  # Skip iteration
            # Case cluster is not in valid ones
            if cluster_name not in cluster_names:
                continue
            # Move cluster directory to target directory
            shutil.move(cluster_path, target_path)

        # Update timer
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('Batch {:s} done,'.format(batch_name), end=' ')
            print('with size {:d}'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Update log
        log({
            'end_time': time_tot,
            'end_size': len(cluster_names)
        })

        # Return tot time and cluster names
        return time_tot, cluster_names

    # Run build
    def make_build(self, build_path, author_name, verbose=False):
        """Make build of HMMs

        First, makes automatic SEED trimming among clusters in given build
        directory. Invalid trimmed clusters will be moved to NOALN/ directory.

        Then, makes HMMs and seaches them against Uniprot using pfbuild, which
        intrinsically calls pfmake either.

        New HMM searches results are then checked to spot UniProt intersections
        and eventually intersectiing clusters are moved from build directory
        to Uniprot/ subfolder.

        Afterwards, new HMMs are build and new MGnifam entries are created by
        running mgnifam (HMM search against MGnifam dataset).

        Lastly, new automatic annotation is made, feeding author name to
        annonation script.

        Args
        build_path (str)        Path to build directory, where clusters folders
                                are actually stored
        author_name (str)       Name of the author which ran the pipeline
        verbose (bool)          Wether to show verbose output or not

        Return
        (float)                 Time taken to run
        (set)                   Set of cluster names which ran successfully
        """

        # Start timer
        time_beg, time_end = time(), 0.0

        # Define iterator through clusters
        clusters_iter = iglob(os.path.join(build_path, 'MGYP*'))
        # Initialize set of cluster names
        cluster_names = set()
        # Loop through each cluster
        for cluster_path in clusters_iter:
            # Get cluster name from path
            cluster_name = os.path.basename(cluster_path)
            # Save cluster name
            cluster_names.add(cluster_name)

        # Verbose log
        if verbose:
            print('Making build of size {:b}'.format(len(cluster_names)))

        # Define NOALN path
        noaln_path = os.path.join(build_path, 'NOALN')
        # Check if NOALN folder exists
        if not os.path.exists(noaln_path):
            # Try to make new NOALN folder
            os.mkdir(noaln_path)

        # Define DISCARD path
        discard_path = os.path.join(build_path, 'DISCARD')
        # Check if DISCARDED path exist
        if not os.path.exists(discard_path):
            # Try to make new DISCARD folder
            os.mkdir(discard_path)

        # Define log path
        log_path = os.path.join(build_path, 'build.json')
        # Initialize log
        log = Log(log_path=log_path, log_dict={
            'build_path': build_path,
            'build_size': len(cluster_names)
        })

        # Make SEED automatic trimming
        time_run, cluster_names = self.trim_seed_alignments(
            clusters_path=build_path,
            noaln_path=noaln_path,
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Update log
        log({
            'trim_time': time_run,
            'trim_size': len(cluster_names)
        })

        # Make HMMs and search them against UniProt
        time_run, cluster_names = self.make_hmm_models(
            clusters_path=build_path,
            cluster_names=cluster_names,
            db='uniprot',
            withpfmake=True,
            make_eval=None,
            verbose=verbose
        )

        # Update log
        log({
            'pfbuild_time': time_run,
            'pfbuild_size': len(cluster_names)
        })

        # Check HMMs for UniProt insersection
        time_run, cluster_names = self.check_hmm_uniprot(
            clusters_path=build_path,
            discard_path=discard_path,
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Update log
        log({
            'uniprot_check_time': time_run,
            'uniprot_check_size': len(cluster_names)
        })

        # Make HMMs and search them against MGnifam
        time_run, cluster_names = self.make_hmm_models(
            clusters_path=build_path,
            cluster_names=cluster_names,
            db='mgnifam',
            withpfmake=True,
            make_eval=0.01,
            verbose=verbose
        )

        # Update log
        log({
            'mfbuild_time': time_run,
            'mfbuild_size': len(cluster_names)
        })

        # Make annotations
        time_run, cluster_names = self.make_annotations(
            clusters_path=build_path,
            discard_path=discard_path,
            cluster_names=cluster_names,
            author_name=author_name,
            verbose=verbose
        )

        # Update log
        log({
            'annotation_time': time_run,
            'annotation_size': len(cluster_names)
        })

        # Update timer
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('Build done,', end=' ')
            print('with size {:d}'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Update log
        log({
            'end_time': time_tot,
            'end_size': len(cluster_names)
        })

        # Return tot time and cluster names
        return time_tot, cluster_names

    # Make SEED alignments
    def make_seed_alignments(self, clusters_path, cluster_names, verbose=False):
        """Make SEED alignments

        Given some cluster names and an output directory path, retrieves SEED
        alignments for those clusters and puts the related folders in output
        folder. Returns set of clustr names which successfully ran.

        Args
        clusters_path (str)     Path to output folder, where to store clusters
                                related folders.
        cluster_names (set)     Set of cluster names (unique) to retrieve from
                                MGnifam dataset
        verbose (bool)          Wether to return boolean output

        Return
        (float)                 Time taken to run
        (set)                   Set of clusters which ran successfully
        """
        # Verbose log
        if verbose:
            print('Running MGSEED for {:d} clusters'. format(len(cluster_names)))

        # Initialize timers
        time_beg, time_end = time(), 0.0

        # Initialize job dict(cluster name: job object)
        bjobs = dict()
        # Loop through each cluster name
        for cluster_name in cluster_names:
            # Run mgseed script
            bjobs[cluster_name] = self.run_mgseed(cluster_name, cwd=clusters_path)
        # Check jobs
        Bjob.check(bjobs.values(), delay=30, verbose=verbose)

        # Get set of cluster names whose status is DONE
        cluster_names = {c for c in bjobs if bjobs.get(c).is_done()}

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('SEED alignments done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Return time taken
        return time_tot, cluster_names

    # Make SEED automatic trimming
    def trim_seed_alignments(self, clusters_path, noaln_path, cluster_names, verbose=False):
        """Automatically trim SEED alignments

        For given cluster names, check that the associated cluster directory
        exists, then takes the old SEED alignment from it and makes a new
        automatically trimmed alignment.

        Retrieved alignment could be empty: in that case, cluster will be
        discarded by moving its folder to NOALN path.

        Args
        clusters_path (str)     Path where clusters folders are stored
        noaln_path (str)        Path where discarded clusters folders will be
                                moved (discarded)
        cluster_names (list)    Names of clusters whose SEED must be trimmed

        Return
        (float)                 Total time taken for script to run
        (set)                   Set of clusters which ran successfully
        """
        # Verbose log
        if verbose:
            print('Trimming SEED alignments', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)))

        # Define timers
        time_beg, time_end = time(), 0.0

        # Check that NOALN directory exist
        if not os.path.exists(noaln_path):
            # Try to make one
            os.mkdir(noaln_path)

        # Define clusters iterator
        clusters_iter = iglob(os.path.join(clusters_path, 'MGYP*'))

        # Loop through each cluster path
        for cluster_path in clusters_iter:
            # Get cluster name out of cluster path
            cluster_name = os.path.basename(cluster_path)
            # Ensure cluster name is a valid one
            if not cluster_name in cluster_names:
                continue
            # Ensure cluster path exists
            if not os.path.exists(cluster_path):
                # Update cluster names set
                cluster_names.remove(cluster_name)
                continue
            # Trim cluster
            self.trim_seed(cluster_path, noaln_path, self.seed_transform, verbose=verbose)

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('Trimming SEED alignments done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Return either time taken and set of cluster names
        return time_tot, cluster_names

    # Make HMMs
    def make_hmm_models(self, clusters_path, cluster_names, verbose=False, **kwargs):
        """Make HMM models

        Build HMM models out of clusters whose name is given as input and which
        have an associated directory in given clusters directory, specified by
        input path.

        Built HMMs are searched against a target dataset by pfmake, a script
        being called from within pfbuild script.

        This method returns only names of clusters for which an HMM has
        been successfully built.

        Args
        clusters_path (str)     Path where cluster folders can be found
        cluster_names (list)    List of cluster names to run through pfbuild
        verbose (bool)          Wether to show verbose log or not
        **kwargs                Other keyword arguments which can be passed to
                                pfbuild script

        Return
        (float)                 Total time taken to run
        (set)                   Names of the clusters which successfully run
        """
        # Verbose log
        if verbose:
            print('Running PFBUILD for {:d} clusters'.format(len(cluster_names)))

        # Define timers
        time_beg, time_end = time(), 0.0

        # Define clusters iterator
        clusters_iter = iglob(os.path.join(clusters_path, 'MGYP*'))

        # Initialize jobs dict(cluster name: bjob object)
        bjobs = dict()
        # Loop thriugh each cluster
        for cluster_path in clusters_iter:
            # Define cluster name
            cluster_name = os.path.basename(cluster_path)
            # Ensure that cluster path exists
            if not os.path.exists(cluster_path):
                continue  # Skip iteration
            # Ensure that SEEED alignment has been succesfully run
            if cluster_name not in cluster_names:
                continue  # Skip iteration
            # Run PFBUILD
            bjobs[cluster_name] = self.run_pfbuild(cluster_path, **kwargs)
        # Check jobs
        Bjob.check(bjobs.values(), delay=30, verbose=verbose)

        # Get set of cluster names whose status is DONE
        cluster_names = {c for c in bjobs if bjobs[c].get_status() == 'DONE'}

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('PFBUILD done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Return time taken and successful clusters names
        return time_tot, cluster_names

    # Check HMMs against uniprot
    def check_hmm_uniprot(self, clusters_path, discard_path, cluster_names, verbose=False):
        """Check HMM results against UniProt

        Checks all HMM searches against UniProt results, discards directories
        associated to clusters with at least one hit in UniProt, by moving them
        into UniProt/ folder.

        Before running the check, this method ensures that only clusters in
        given set of cluster names will be actually run, hence moves any other
        cluster into DISCARDED folder, specified by given discard path.

        Args
        clusters_path (str)     Path where cluster folders are stored
        discard_path (str)      Path where clusters not willing to be checked
                                must be moved to
        cluster_names (set)     Set of clusters which must be run
        verbose (bool)          Whether to show verbose log or not
        """
        # Verbose log
        if verbose:
            print('Checking {:d}'.format(len(cluster_names)), end=' ')
            print('HMM models against UniProt')

        # Initialize timers
        time_beg, time_end = time(), 0.0

        # Initialize cluster iterator
        clusters_iter = iglob(os.path.join(clusters_path, 'MGYP*'))

        # Discard clusters
        discarded_names, cluster_names = self.discard_clusters(
            clusters_path=clusters_path,
            discard_path=discard_path,
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Run check uniprot
        self.check_uniprot(clusters_path)

        # Reinitialize cluster names
        cluster_names = set()
        # Go through clusters again
        for cluster_path in clusters_iter:
            # Define cluster name
            cluster_name = os.path.basename(cluster_path)
            # Ensure cluster path still exist (not moved to uniprot)
            if not os.path.exists(cluster_path):
                continue
            # Add cluster name to set of valid ones
            cluster_names.add(cluster_name)

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('HMM checked against Uniprot', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Return time taken
        return time_tot, cluster_names

    def make_annotations(self, clusters_path, discard_path, cluster_names, author_name, verbose=False):
        """Annotate clusters

        Discard clusters not in given set of cluster names, then makes
        annotations for the remaining ones.

        Args
        clusters_path (str)     Path to directory where clusters are stored
        discard_path (str)      Path where to move not valid clusters
        cluster_names (set)     Set of valid cluster names
        author_name (str)       Name of the annotator
        verbose (bool)          Wether to show verbose log or not

        Return
        (float)                 Time taken to run
        (set)                   Set of remaining cluster names
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0

        # Verbose log
        if verbose:
            print('Annotating {:d} clusters'.format(len(cluster_names)))

        # Initialize cluster iterator
        clusters_iter = iglob(os.path.join(clusters_path, 'MGYP*'))

        # Define names of discarded clusters
        discarded_names, cluster_names = self.discard_clusters(
            clusters_path=clusters_path,
            discard_path=discard_path,
            verbose=verbose
        )

        # Run annotation script
        self.run_annotate()

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('Annotated {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.2f} seconds'.format(time_tot))

        # Return time taken to run and set of valid cluster names
        return time_tot, cluster_names

    # Discard clusters whse name is not in given set
    def discard_clusters(self, clusters_path, discard_path, cluster_names, verbose=False):
        """Discard clusters

        Takes as input a set of cluster names, then discards those cluster
        folders whose name is not in given set, by moving them from their
        actual path to given discard path.

        Args

        verbose (bool)      Wether to show verbose output or not

        Return
        (set)               Set of cluster names discarded
        (set)               Set of cluster names not discarded
        """
        # Initialize cluster iterator
        clusters_iter = iglob(os.path.join(clusters_path, 'MGYP*'))

        # Define names of discarded clusters
        discarded_names = set()
        # Discard clusters
        for cluster_path in clusters_iter:
            # Define cluster name from its path
            cluster_name = os.path.basename(cluster_path)
            # Ensure cluster name is not valid
            if cluster_name in cluster_names:
                continue
            # Move not valid cluster to discard path
            shutil.move(cluster_path, os.path.join(discard_path, cluster_name))
            # Add name to discarded ones
            discarded_names.add(cluster_name)
            # Remove cluster name from valid ones
            cluster_names.remove(cluster_name)

        # Verbose log
        if verbose:
            print('Discarded {:d} clusters, kept {:d} clusters'.format(
                len(discarded_names), len(cluster_names)
            ))

        # Return set of cluster names
        return discarded_names, cluster_names

    # Run mgseed.pl
    def run_mgseed(self, cluster_name, cwd='./'):
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
    def run_pfbuild(self, cluster_path, db='uniprot', withpfmake=True, make_eval=None):
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

    # # Run pfbuild with mfbuild parameters as default
    # def mfbuild(self, cluster_path, db='mgnify', withpfmake=True, make_eval=0.01):
    #     # Just pass parameters to pfbuild script
    #     return self.pfbuild(cluster_path=cluster_path, db=db, withpfmake=withpfmake, make_eval=make_eval)

    # Run check_uniprot.pl
    def check_uniprot(self, clusters_path):
        # Run check_uniprot.pl in current batch directory
        return subprocess.run(
            capture_output=True,
            encoding='utf-8',
            check=True,
            env=self.env,
            cwd=clusters_path,
            args=['check_uniprot.pl']
        )

    # Annotate MGnifam release
    def run_annotate(self, clusters_path, author_name):
        return subprocess.run(
            capture_output=True,
            encoding='utf-8',
            check=True,
            env=self.env,
            cwd=clusters_path,
            args=['annotateMGDUF.pl', '-directory', '.', '-author', author_name]
        )

    # Automatically trim SEED alignment
    def trim_seed(self, cluster_path, noaln_path, transform, verbose=False):
        """Trim alignments and plot results

        Duplicates SEED alignment in given cluster path, then trims one making
        it the sew SEED alignment, while the other will be stored as raw one
        but never used in later computations.

        In case automatic trimming of SEED alignment produces an empty result,
        cluster will be moved to given NOALN path, thresefore discarding it.

        Args
        cluster_path (str)      Path to cluster whose SEED alignment must be
                                trimmed
        noaln_path (str)        Path where empty trimmed alignments will be
                                moved (discarded)
        transform (Transform)   Transformer which will be fed with original
                                SEED alignment to produce the trimmed one
        verbose (bool)          Wether to show verbose log or not
        """
        # Get cluster name from cluster path
        cluster_name = os.path.basename(cluster_path)

        # Define discard path
        discard_path = os.path.join(noaln_path, cluster_name)

        # Define SEED paths
        raw_path = os.path.join(cluster_path, 'SEED_raw')
        new_path = os.path.join(cluster_path, 'SEED')

        # Rename original cluster SEED alignment
        shutil.move(new_path, raw_path)

        # Load raw SEED multiple sequence alignment from file
        raw_seed = MSA.from_aln(raw_path)
        # Define seed shape
        raw_shape = raw_seed.aln.shape
        # Verbose log
        if verbose:
            print('Raw seed alignment has shape {:d} x {:d}'.format(*raw_shape))

        # Go through multiple sequence alignment transformation pipeline
        new_seed = transform(raw_seed)

        # Case trimmed alignment is empty
        if new_seed.is_empty():
            # Discard current cluster
            shutil.move(cluster_path, discard_path)
            # Verbose log
            if verbose:
                print('Cluster {:s} was empty, then has been discarded'.format(
                    cluster_name
                ))

        # Otherwise
        else:
            # Store new SEED multiple sequence alignment to file
            new_seed.to_aln(new_path)
            # Get new shape
            new_shape = new_seed.aln.shape
            # Verbose log
            if verbose:
                print('New SEED alignment for cluster {:s}'.format(cluster_name), end=' ')
                print('has shape {:d} x {:d}'.format(*new_shape))


# Main
if __name__ == '__main__':

    # Setup local environment
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

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Make new MGnifam release'
    )
    parser.add_argument(
        '--in_path', nargs='+', type=str, required=True,
        help='Path to one or more files containing cluster names'
    )
    parser.add_argument(
        '--out_path', type=str, required=True,
        help='Output directory path of the release'
    )
    parser.add_argument(
        '-a', '--author_name', type=str, required=True,
        help='Author name which is making the new release'
    )
    parser.add_argument(
        '-n', '--max_clusters', type=int, default=1000,
        help='Maximum number of clusters to run through the pipeline'
    )
    parser.add_argument(
        '-b', '--batch_size', type=int, default=100,
        help='Maximum number of clusters per batch to run'
    )
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Whether to show verbose log or not'
    )
    args = parser.parse_args()

    # Instantiate new pipeline
    pipeline = OldRelease(delay=30, env=env)

    # Run pipeline
    pipeline(
        in_path=args.in_path,
        out_path=args.out_path,
        author_name=args.author_name,
        batch_size=args.batch_size,
        max_clusters=args.max_clusters,
        verbose=bool(args.verbose)
    )
