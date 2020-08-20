# Dependencies
from src.pipeline.pipeline import Pipeline
from src.scheduler import LSFScheduler, LocalScheduler
from src.dataset import LinClust, Fasta
from src.sequences import fasta_iter
from src.disorder import MobiDbLite

from src.msa.transform import Compose, OccupancyTrim, OccupancyFilter
from src.msa.msa import MSA, Muscle
from src.hmm.hmmbuild import HMMBuild
from src.hmm.hmmsearch import HMMSearch
from src.hmm.hmmalign import HMMAlign
from src.hmm.hmmer import Tblout, Domtblout
from src.hmm.hmmer import merge_scores, merge_hits
from src.hmm.hmmer import to_tsv
from src.hmm.hmm import HMM

from src.mgnifam import Cluster
from src.utils import is_gzip, gunzip, as_bytes
from subprocess import CalledProcessError
from tempfile import mkdtemp, mkstemp
from glob import glob
from time import time
import argparse
import shutil
import random
import json
import math
import os
import re


class Build(Pipeline):
    """ Make a batch of LinCLust clusters

    This pipeline takes as input a list of cluster names, for each cluster
    retrieves its member sequences, makes a SEED alignment through Muscle
    multiple sequence alignment algorithm, trims it automatically (and
    eventually discards it, if not valid)
    """

    # Constructor
    def __init__(
        # Pipeline parameters, required to handle job scheduling
        # self, cluster_type, cluster_kwargs,
        self, scheduler,
        # Path to datasets (fixed)
        linclust_path, mgnifam_path, uniprot_path,
        # Compositional bias threshold settings
        comp_bias_threshold=0.2, comp_bias_inclusive=True,
        # Automatic trimming settings
        trim_threshold=0.4, trim_inclusive=True,
        filter_threshold=0.5, filter_inclusive=True,
        # Post trimming settings
        seed_min_width=1, seed_min_height=1,
        # Search against UniProt settings
        uniprot_e_value=0.01, uniprot_height=None, uniprot_width=None,
        # Search against MGnifam settings
        mgnifam_e_value=0.01, mgnifam_height=None, mgnifam_width=None,
        # Command line arguments
        mobidb_cmd=['mobidb_lite.py'],  # Path to MobiDB Lite predictor
        muscle_cmd=['muscle'],  # Path to Muscle alignment algorithm
        hmmbuild_cmd=['hmmbuild'],  # Path to hmmbuild script
        hmmalign_cmd=['hmmalign'],  # Path to hmmalign script
        hmmsearch_cmd=['hmmsearch'],  # Path to hmmsearch script
        # Environmental variables
        env=os.environ.copy()
    ):
        # Save datasets path
        self.linclust = LinClust.from_list(linclust_path)
        self.mgnifam = Fasta.from_list(mgnifam_path)
        self.uniprot = Fasta.from_list(uniprot_path)
        # Save compositional bias parameters
        self.comp_bias_threshold = comp_bias_threshold
        self.comp_bias_inclusive = comp_bias_inclusive
        # Save automatic trimming settings
        self.trim_threshold = trim_threshold
        self.trim_inclusive = trim_inclusive
        self.filter_threshold = filter_threshold
        self.filter_inclusive = filter_inclusive
        # Save post trimming settings
        self.seed_min_width = seed_min_width
        self.seed_min_height = seed_min_height
        # Define automatic trimming
        self.seed_transform = Compose([
            # Exclude regions outside N- and C- terminal
            OccupancyTrim(
                threshold=trim_threshold,
                inclusive=trim_inclusive
            ),
            # Exclude sequences with less than half occupancy
            OccupancyFilter(
                threshold=filter_threshold,
                inclusive=filter_inclusive
            )
        ])
        # Save UniProt hmmsearch settings
        self.uniprot_e_value = uniprot_e_value
        self.uniprot_height = uniprot_height
        self.uniprot_width = uniprot_width
        # Save MGnifam hmmsearch settings
        self.mgnifam_e_value = mgnifam_e_value
        self.mgnifam_height = mgnifam_height
        self.mgnifam_width = mgnifam_width
        # Define scripts interfaces
        self.mobidb = MobiDbLite(cmd=mobidb_cmd, env=env)
        self.muscle = Muscle(cmd=muscle_cmd, env=env)
        self.hmmbuild = HMMBuild(cmd=hmmbuild_cmd, env=env)
        self.hmmalign = HMMAlign(cmd=hmmalign_cmd, env=env)
        self.hmmsearch = HMMSearch(cmd=hmmsearch_cmd, env=env)
        # Save scheduler
        self.scheduler = scheduler
        # Store environmental variables
        self.env = env

    # Run the pipeline
    def run(self, cluster_names, clusters_path, author_name='', min_jobs=1, max_jobs=100, verbose=False):
        """ Run BATCH pipeline

        Args
        cluster_names (set)         Set of cluster names whose SEED alignment
                                    must be generated
        clusters_path (str)         Path where to store results, an attempt of
                                    making the output directory is made if it
                                    does not exist yet.
        min_jobs (int)              Minimum number of jobs to spawn
        max_jobs (int)              Maximum number of jobs to spawn
        verbose (bool)              Whether to print verbose log or not

        Return

        Raise
        """

        # Verbose log
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Running build pipeline', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('at {:s}'.format(clusters_path))

        # Initialize sub-directories
        bias_path, noaln_path, uniprot_path = self.init_subdir(clusters_path)

        # Adapt scheduler
        self.scheduler.adapt(min_jobs, max_jobs)
        # Open new client
        with self.scheduler:

            # Fill cluster members dictionary
            cluster_members = self.search_cluster_members(
                cluster_names=cluster_names,
                client=self.scheduler.client,
                verbose=verbose
            )

            # Check cluster members
            self.check_cluster_members(cluster_names, cluster_members)

            # Initialize set of retrieved sequence accession numbers
            sequences_acc = set()
            # Loop through each cluster
            for members_acc in cluster_members.values():
                # Update sequence accession numbers
                sequences_acc |= set(members_acc)

            # Fill mgnifam dictionary
            fasta_sequences, _ = self.search_fasta_sequences(
                sequences_acc=sequences_acc,
                fasta_dataset=self.mgnifam,
                client=self.scheduler.client,
                verbose=verbose
            )

            # Check fasta sequences
            self.check_fasta_sequences(sequences_acc, fasta_sequences)

            # Make clusters directories and SEED.fa fasta files
            self.make_seed_sources(
                clusters_path=clusters_path,
                cluster_members=cluster_members,
                fasta_sequences=fasta_sequences,
                verbose=verbose
            )

            # Check compositional bias
            self.check_comp_bias(
                clusters_path=clusters_path,
                bias_path=bias_path,
                client=self.scheduler.client,
                verbose=verbose
            )

            # Make SEED alignments
            self.make_seed_alignments(
                clusters_path=clusters_path,
                noaln_path=noaln_path,
                client=self.scheduler.client,
                verbose=verbose
            )

            # Make HMM models
            self.make_hmm_models(
                clusters_path=clusters_path,
                verbose=verbose,
                client=self.scheduler.client
            )

            # Check if UniProt shape hasn't been already set
            if not (self.uniprot_height and self.uniprot_width):
                # Retrieve UniProt size and it longest sequence length
                _, self.uniprot_width, self.uniprot_height = self.get_longest(
                    dataset_chunks=self.uniprot,
                    verbose=verbose,
                    client=self.scheduler.client
                )

            # Retrieve UniProt shape
            uniprot_shape = (self.uniprot_height, self.uniprot_width)
            # Verbose
            if verbose:
                print('UniProt dataset has shape', end=' ')
                print('{:d} x {:d}'.format(*uniprot_shape))

            # Check if MGnifam shape hasn't been already set
            if not (self.mgnifam_height and self.mgnifam_width):
                # Retrieve MGnifam size and it longest sequence length
                _, self.mgnifam_width, self.mgnifam_height = self.get_longest(
                    dataset_chunks=self.mgnifam,
                    verbose=verbose,
                    client=self.scheduler.client
                )

            # Retrieve MGnifam shape
            mgnifam_shape = (self.mgnifam_height, self.mgnifam_width)
            # Verbose log
            if verbose:
                print('MGnifam dataset has shape', end=' ')
                print('{:d} x {:d}'.format(*mgnifam_shape))

        # Search against UniProt pipeline
        self.search_hmm_uniprot(
            clusters_path=clusters_path,
            uniprot_path=uniprot_path,
            uniprot_shape=uniprot_shape,
            min_jobs=min_jobs,
            max_jobs=max_jobs,
            verbose=verbose
        )

        # Run search against MGnifam sub-pipeline
        self.search_hmm_mgnifam(
            clusters_path=clusters_path,
            mgnifam_shape=mgnifam_shape,
            author_name=author_name,
            min_jobs=min_jobs,
            max_jobs=max_jobs,
            verbose=verbose
        )

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('Build pipeline done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.0f} seconds'.format(time_end - time_beg))

        # Exit
        return

    # Initialize sub-directories
    def init_subdir(self, clusters_path):
        """ Initialize sub-directories

        Makes three main sub folders, if not already there: bias (for clusters
        discarded due to high bias), noaln (for clusters discarded due to
        invalid alignment), discard (for cluster discarded for other reasons)
        and uniprot (for clusters discarded due to UniProt intersection)

        Args
        clusters_path (str)     Path to clusters root directory

        Return
        (str)                   Path to bias/ subdirectory
        (str)                   Path to noaln/ subdirectory
        (str)                   Path to uniprot/ subdirectory
        """
        # Initialize list of sub-directories
        sub_paths = ['BIAS', 'NOALN', 'UNIPROT']
        # Loop through each sub-directory path
        for i, sub_path in enumerate(sub_paths):
            # Make full sub-directory path
            sub_path = os.path.join(clusters_path, sub_path)
            # Case directory does not already exist
            if not os.path.exists(sub_path):
                # Attempt to make one
                os.mkdir(sub_path)
            # Substitute partial with full sub-directory path
            sub_paths[i] = sub_path
        # Return sub-directories paths
        return tuple(sub_paths)

    # Retrieve longest sequence of a sequences dataset
    def get_longest(self, dataset_chunks, client, verbose=False):
        """ Retrieve longest sequence and dataset length

        Args
        dataset_chunks (list)       List of target dataset chunks to go through
        client (dask.Client)        Dask client used to distribute search
        verbose (bool)              Whether to print out verbose log

        Return
        (str)                       Longest sequence in FASTA dataset
        (int)                       Length of longest sequence in FASTA dataset
        (int)                       Dataset length
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show process start
            print('Searching for longest sequence against', end=' ')
            print('{:d} chunks...'.format(len(dataset_chunks)), end=' ')

        # Initailize list of futures
        futures = list()
        # Loop through each chunk in dataset
        for chunk in dataset_chunks:
            # Submit length function
            future = client.submit(chunk.get_longest)
            # Store future
            futures.append(future)

        # Retrieve results
        results = client.gather(futures)

        # Initialize longest sequence and its length
        longest_seq, longest_len = '', 0
        # Initialize dataset length
        dataset_size = 0
        # Loop through each result
        for i in range(len(results)):
            # Get either longest sequence and dataset length
            chunk_longest_seq, chunk_longest_len, chunk_size = results[i]
            # Case current sequence is longer than previous longest
            if longest_len < chunk_longest_len:
                # Update longest sequences
                longest_seq = chunk_longest_seq
                longest_len = chunk_longest_len
            # Update dataset length
            dataset_size += chunk_size

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Return both longest sequence and dataset length
        return longest_seq, longest_len, dataset_size

    # Search for cluster members
    def search_cluster_members(self, cluster_names, client, verbose=False):
        """ Search for cluster member sequences

        Args
        cluster_names (set)         Set of cluster names to search against
                                    LinClust dataset
        client (Client)             Dask Client used to schedule jobs
        verbose (bool)              Whether to print verbose log or not

        Return
        (dict)                      Mapping cluster names (keys) to set of
                                    sequence accession numbers (values)
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Verbose
        if verbose:
            print('Searching member sequences for', end=' ')
            print('{:d} clusters...'.format(len(cluster_names)), end=' ')

        # Initialize futures
        futures = list()
        # Go through each LinClust dataset chunk
        for chunk in self.linclust:
            # Submit search function to client
            future = client.submit(chunk.search, cluster_names=cluster_names)
            # Store future
            futures.append(future)

        # Retrieve distributed computation results
        results = client.gather(futures)

        # Initialize cluster members dictionary (output)
        cluster_members = dict()
        # Loop through resulting clusters dictionaries
        for i in range(len(results)):
            # Loop through every cluster name in result dictionary
            for cluster_name in results[i].keys():
                # Initialize entry for current cluster
                cluster_members.setdefault(cluster_name, set())
                # Update cluster members list
                cluster_members[cluster_name] |= set(results[i][cluster_name])

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            # Print execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Return cluster members
        return cluster_members

    # Check cluster members found
    def check_cluster_members(self, cluster_names, cluster_members):
        # Check that cluster members keys are the same as cluster names
        if set(cluster_names) != set(cluster_members.keys()):
            # Raise new exception and interupt execution
            raise KeyError(' '.join([
                'Error: cluster members have not been retrieved for',
                'all the {:d} input clusters: '.format(len(cluster_names)),
                'aborting pipeline execution'
            ]))

    # Search for fasta sequences
    def search_fasta_sequences(self, sequences_acc, fasta_dataset, client, verbose=False):
        """ Retrieve fasra sequences from a Fasta dataset

        Args
        sequences_acc (set)         Set of sequence accession numbers to search
                                    for against given target FASTA dataset
        fasta_dataset (Dataset)     Target FASTA dataset instance, whose search
                                    method will be used in order to find
                                    FASTA entries associated to given
                                    accession numbers
        client (Client)             Dask Client used to schedule jobs
        verbose (bool)              Whether to print verbose log or not

        Return
        (dict)                      Dictionary mapping sequence accession
                                    numbers (keys) to list of fasta entries
                                    (values)
        (int)                       Total length of input FASTA dataset
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Verbose log
        if verbose:
            print('Searching {:d}'.format(len(sequences_acc)), end=' ')
            print('fasta sequences...', end=' ')

        # Initialize futures
        futures = list()
        # Loop through each chunk in given FASTA dataset
        for chunk in fasta_dataset:
            # Submit search function against FASTA dataset
            future = client.submit(chunk.search, sequences_acc, ret_length=True)
            # Save future
            futures.append(future)

        # Gather results
        results = client.gather(futures)

        # Initialize FASTA sequences dict(sequence accession: fasta entry)
        fasta_sequences = dict()
        # Initialize length of FASTA file
        fasta_length = 0
        # Go through each result
        for i in range(len(results)):
            # Retrieve either chunk sequences dict and chunk length
            chunk_sequences, chunk_length = results[i]
            # Update cluster sequences dictionary
            for sequence_acc in chunk_sequences:
                fasta_sequences[sequence_acc] = chunk_sequences[sequence_acc]
            # Update length of FASTA file
            fasta_length = fasta_length + chunk_length

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            # Print time
            print('done in {:.0f} seconds'.format(time_end - time_beg))
            # Print results size
            print('Retrieved {:d} fasta'.format(len(fasta_sequences)), end=' ')
            print('entries from dataset of size {:d}'.format(fasta_length))

        # Return either the sequences dictionary and the dataset length
        return fasta_sequences, fasta_length

    # Check fasta sequences found
    def check_fasta_sequences(self, sequences_acc, fasta_sequences):
        # Loop through all required seuqnece accession numbers
        for sequence_acc in sequences_acc:
            # Check if sequnece accession is available and not empty
            if not fasta_sequences.get(sequence_acc, ''):
                # Raise new error
                raise KeyError(' '.join([
                    'Error: required sequence accession number',
                    '{} has no associated fasta entry:'.format(sequence_acc),
                    'aborting pipeline execution'
                ]))

    # Make clusters directories and write SEED.fa files
    def make_seed_sources(self, clusters_path, cluster_members, fasta_sequences, verbose=False):
        """ Make cluster source files (SEED.fa)

        Loop through all cluster names (keys) in cluster members dictionary:
        for each cluster member, retrieve the associated fasta entry in given
        fasta dictionary.

        If a cluster directory does not exist, tries to make a new one.

        Args
        clusters_path (str)         Path where clusters directory can be found
        cluster_members (dict)      Dictionary mapping cluster names to cluster
                                    members (i.e. sequence accession numbers)
        fasta_sequences (dict)      Dictionary mapping sequence accession
                                    numbers to fasta entries
        verbose (bool)              Whether to show verbose log

        Raise
        (OSError)                   In case problems occurred while making
                                    new directories and files
        """
        # Loop through each cluster name
        for cluster_name in cluster_members:
            # Try to instantiate directory
            try:
                # Define cluster path
                cluster_path = os.path.join(clusters_path, cluster_name)
                # Define SEED.fa file path
                fasta_path = os.path.join(cluster_path, 'SEED.fa')

                # Check if directory already exists
                if not os.path.isdir(cluster_path):
                    # Try making new directory
                    os.mkdir(cluster_path)

                # Make SEED.fa file
                with open(fasta_path, 'w') as fasta_file:
                    # Loop through each sequence accession in clurrent cluster
                    for acc in cluster_members[cluster_name]:
                        # Retrieve fasta entry from sequences dictionary
                        fasta_entry = fasta_sequences[acc]
                        # Write fasta entry to file
                        fasta_file.write(fasta_entry + '\n')

            except OSError:
                # Raise
                raise OSError(' '.join([
                    'could not make cluster at',
                    '{:s}:'.format(cluster_path),
                    'aborting execution'
                ]))

        # Verbose
        if verbose:
            print('Source SEED.fa done for all the', end=' ')
            print('{:d} clusters'.format(len(cluster_members)))

    # Check compositional bias
    def check_comp_bias(self, clusters_path, bias_path, client, verbose=False):
        """ Check compositional bias

        Loop through all clusters in given clusters path and runs MobiDB Lite
        predictor over SEED.fa file (if any). Then computes compositional
        bias over prediction results: if too high, then discard the cluster
        sub-directory by moving it to bias folder.

        Disorder predictions will be stored inside the cluster sub-folder
        in mobidb.disorder file, while compositional bias value will be stored
        in mobidb.bias file

        Args
        clusters_path (str)         Path where clusters directory can be found
        client (Client)             Dask client used to submit jobs
        verbose (bool)              Wether to print verbose output or not

        Return
        (dict)                      Dictionary mapping cluster name to
                                    compositional bias
        """
        # Intialize timers
        time_beg, time_end = time(), 0.0
        # Define clusters iterator
        clusters_iter = glob(os.path.join(clusters_path, 'MGYP*'))
        # Verbose log
        if verbose:
            print('Computing compositional bias for', end=' ')
            print('{:d} clusters...'.format(len(clusters_iter)), end=' ')

        # Initialize futures dict(cluster name: compositional bias)
        futures = dict()
        # Loop through each cluster
        for cluster_path in clusters_iter:
            # Get cluster name
            cluster_name = os.path.basename(cluster_path)
            # Define fasta file path
            fasta_path = os.path.join(cluster_path, 'SEED.fa')
            # Check if SEED.fa file exists
            if not os.path.isfile(fasta_path):
                # Verbose
                if verbose:
                    print('Could not find source SEED.fa file', end=' ')
                    print('for cluster {:s}'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(cluster_path))
                # Skip iteration
                continue

            # Run distributed compositional bias computation
            futures[cluster_name] = client.submit(
                func=self.predict_comp_bias,
                fasta_path=fasta_path,  # Input file
                pred_threshold=1,  # Set threshold for disorder prediction
                pred_inclusive=True,  # Set thresold inclusive
                mobidb=self.mobidb
            )

        # Retrieve results dict(cluster name: comp bias)
        clusters_bias = client.gather(futures)
        # Retrieve compositional bias threshold
        threshold = self.comp_bias_threshold
        inclusive = self.comp_bias_inclusive
        # Loop through each cluster in cluster iterator
        for cluster_path in clusters_iter:
            # Retrieve cluster name
            cluster_name = os.path.basename(cluster_path)
            # Retrieve cluster bias
            cluster_bias = clusters_bias[cluster_name]

            # Case bias threshold is inclusive
            if inclusive and cluster_bias < threshold:
                # Skip iteration
                continue
            # Case bias threshold is exclusive
            if not inclusive and cluster_bias <= threshold:
                # Skip iteration
                continue
            # Move cluster to bias directory
            shutil.move(cluster_path, os.path.join(bias_path, cluster_name))

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

            # Compute number of clusters above that threshold
            clusters_above = sum([
                int((b >= threshold and inclusive) or (b > threshold))
                for b in clusters_bias.values()
            ])
            # Compute number of clusters below that threshold
            clusters_below = sum([
                int((b < threshold and inclusive) or (b <= threshold))
                for b in clusters_bias.values()
            ])
            # # Show compositional bias
            print('Compositional bias threshold is {:.2f}:'.format(threshold))
            print('  {:d} clusters are above it'.format(clusters_above))
            print('  {:d} clusters are below it'.format(clusters_below))

    @staticmethod
    def predict_comp_bias(fasta_path, mobidb, pred_threshold=1, pred_inclusive=True):

        # Define path to current cluster directory
        cluster_path = os.path.dirname(fasta_path)
        # Define path to output file
        out_path = os.path.join(cluster_path, 'MOBIDB.out')
        # Define path to BIAS file
        bias_path = os.path.join(cluster_path, 'MOBIDB.bias')

        # Define fasta sequences dict(acc: fasta entry)
        fasta_sequences = dict()
        # Open fasta file
        with open(fasta_path, 'r') as fasta_file:
            # Iterate fasta sequences
            for entry in fasta_iter(fasta_file):
                # Split fasta entry into header and residues
                head, _ = tuple(entry.split('\n'))
                # Retrieve accession out of header
                acc = re.search(r'^>(\S+)', head).group(1)
                # Store fasta sequence
                fasta_sequences.setdefault(acc, entry)

        # Just run MobiDB Lite predictor
        mobidb.run(fasta_path=fasta_path, out_path=out_path)
        # Parse results
        disorder = mobidb.parse(out_path)
        # Compute compositional bias
        comp_bias = disorder.comp_bias(
            fasta_sequences=fasta_sequences,
            threshold=pred_threshold,
            inclusive=pred_inclusive
        )

        # Write out bias file
        with open(bias_path, 'w') as bias_file:
            # Write out bias
            bias_file.write('{:f}'.format(comp_bias))

        # Return compositional bias
        return comp_bias

    # Make SEED alignments
    def make_seed_alignments(self, clusters_path, noaln_path, client, verbose=False):
        """ Make SEED alignments

        Loop through all cluster directories, retrieve source fasta
        sequences file SEED.fa and use it to make a new alignment.
        Then, automatically trim the alignment and check its dimension: if it
        is too small, drop it.

        Args
        clusters_path (str)         Path where cluster folders can be found
        client (str)                Client used to schedule jobs
        verbose (bool)              Wether to show verbose log or not
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Retrieve all clusters in given root directory
        clusters_iter = glob(os.path.join(clusters_path, 'MGYP*'))

        # Verbose
        if verbose:
            print('Making SEED alignments', end=' ')
            print('for {:d} clusters'.format(len(clusters_iter)), end=' ')
            print('at {:s}...'.format(clusters_path), end='')

        # Initialize futures
        futures = dict()
        # Loop through each cluster
        for cluster_path in clusters_iter:
            # Define current cluster name
            cluster_name = os.path.basename(cluster_path)
            # Define fasta file path (SEED.fa)
            fasta_path = os.path.join(cluster_path, 'SEED.fa')

            # # Define raw alignment path
            # raw_aln_path = os.path.join(cluster_path, 'SEED_raw')
            #
            # # Initialize cluster sequences
            # cluster_sequences = list()
            # # Loop through each sequence accession in cluster
            # for sequence_acc in cluster_members.get(cluster_name):
            #     # Save fasta sequence
            #     cluster_sequences.append(fasta_sequences[sequence_acc])

            # Make SEED alignment in cluster directory
            futures[cluster_name] = client.submit(
                # Call wrapper function for MUSCLE alignment
                func=self.make_seed_alignment,
                # Path to source fasta file
                fasta_path=fasta_path,
                # MSA transformer instance
                transform=self.seed_transform,
                # Muscle MSA aligner
                muscle=self.muscle,
            )

        # Retrieve results
        clusters_shape = client.gather(futures)
        # Loop through all clusters
        for cluster_name in clusters_shape:
            # Define cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Retrieve current alignment shape
            h, w = clusters_shape[cluster_name]

            # Case new SEED does fulfill requirements
            if (h >= self.seed_min_height) and (w >= self.seed_min_width):
                # Go to next iteration
                continue

            # Otherwise, discard cluster
            shutil.move(cluster_path, os.path.join(noaln_path, cluster_name))

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    @staticmethod
    def make_seed_alignment(fasta_path, transform, muscle):
        """ Use MUSCLE to make a multiple sequence alignment

        Args
        fasta_path (str)                Path to file which must be aligned
        transform (msa.Transform)       MSA transformations which aligned MSA
                                        must undergo
        muscle (msa.Muscle)             Instance of MUSCLE interface object

        Return
        (int)                           Height of the alignment
        (int)                           Witdh of the alignment
        """
        # Define current cluster path
        cluster_path = os.path.dirname(fasta_path)
        # Define MUSCLE output alignment path (overwrite input fasta)
        out_path = os.path.join(cluster_path, 'SEED.fa')
        # Define trimmed alignment path
        aln_path = os.path.join(cluster_path, 'SEED.aln')

        # Run MUSCLE, get raw output
        muscle(
            # Input fasta file path
            fasta_path=fasta_path,
            # Output alignment path
            out_path=out_path,
            # Set quiet execution
            quiet=True
        )

        # Load SEED alignment from fasta output file
        seed = MSA.from_fasta(out_path, regex=r'>(\S+)')
        # Transform SEED alignment
        seed = transform(seed)
        # Save SEED alignment to file
        seed.to_aln(aln_path)

        # Return height and width of the alignment
        return seed.shape

    # Make HMM models
    def make_hmm_models(self, clusters_path, client, verbose=False):
        """ Build HMM models

        For each cluster in given clusters path, create a new HMM model through
        hmmbuild script, then store it in the associated cluster directory.

        Args
        clusters_path (str)     Path where cluster directories can be found
        client (Client)         Client used to spawn jobs on cluster
        verbose (bool)          Whether to print verbose log or not
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Define SEED.aln files
        seed_iter = glob(os.path.join(clusters_path, 'MGYP*', 'SEED.aln'))

        # Verbose
        if verbose:
            print('Making HMM models for', end=' ')
            print('{:d} clusters...'.format(len(seed_iter)), end=' ')

        # Initialize futures
        futures = dict()
        # Loop through each given cluster name
        for seed_path in seed_iter:
            # Define cluster path
            cluster_path = os.path.dirname(seed_path)
            # Define HMM path
            hmm_path = os.path.join(cluster_path, 'HMM.model')

            # # Ensure that cluster directory exists
            # if not os.path.exists(cluster_path):
            #     # Verbose log
            #     if verbose:
            #         print('Cluster {:s} not found'.format(cluster_name), end=' ')
            #         print('at {:s}: skipping it'.format(cluster_path))
            #     # Skip iteration
            #     continue

            # # Ensure that cluster SEED alignment exists
            # if not os.path.exists(seed_path):
            #     # Verbose log
            #     if verbose:
            #         print('SEED alignment for cluster', end=' ')
            #         print('{:s} not found'.format(cluster_name), end=' ')
            #         print('at {:s}: skipping it'.format(cluster_path))
            #     # Skip iteration
            #     continue

            # Define cluster name
            cluster_name = os.path.basename(cluster_path)
            # Submit new job making HMM model
            futures[cluster_name] = client.submit(
                func=self.hmmbuild.run,
                msa_path=seed_path,
                hmm_path=hmm_path,
                name=cluster_name
            )

        # Wait for results
        client.gather(futures)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Print execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    # Make HMM library containing all HMM models
    def make_hmm_library(self, clusters_path, verbose=False):
        """ Concatenate every HMM model in one HMM library

        Args
        clusters_path (str)     Path where cluster directories can be found
        verbose (bool)          Whether to print verbose log or not

        Return
        (str)                   Path to HMM library
        (int)                   Number of HMM models wrapper in HMM library
        (int)                   Maximum HMM size in HMM library
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Define HMM models iterator
        models_iter = glob(os.path.join(clusters_path, 'MGYP*', 'HMM.model'))
        # Define HMM library path
        library_path = os.path.join(clusters_path, 'HMM.library')

        # Verbose log
        if verbose:
            print('Making HMM library out', end=' ')
            print('of {:d} HMM models'.format(len(models_iter)), end=' ')
            print('at {:s}'.format(library_path), end=' ')

        # Define total number of models in library
        library_len = len(models_iter)
        # Intialize maximum model size
        library_width = 0
        # Open output file
        with open(library_path, 'w') as library_file:
            # Loop through each HMM model file
            for model_path in models_iter:
                # Load model, retireve its length
                model = HMM.from_file(model_path)
                # Store length if maximum
                library_width = max(library_width, model.length)
                # Open input file
                with open(model_path, 'r') as model_file:
                    # Copy model to library
                    shutil.copyfileobj(model_file, library_file)

        # # Initialize output HMM library path
        # hmmlib_path = os.path.join(clusters_path, 'HMMLIB')
        # # Open file in write mode
        # hmmlib_file = open(hmmlib_path, 'w')
        #
        # # Initialize names of clusters loaded into HMM library
        # in_hmmlib = set()
        # # Initialize max HMM model model length in MM library
        # hmmlib_len = 0
        # # Loop through each cluster name in given cluster names set
        # for cluster_name in cluster_names:
        #
        #     # Define cluster path
        #     cluster_path = os.path.join(clusters_path, cluster_name)
        #     # Check if cluster path exists
        #     if not os.path.isdir(cluster_path):
        #         # Verbose log
        #         if verbose:
        #             print('Could not find subfolder', end=' ')
        #             print('for cluster {:s}'.format(cluster_name), end=' ')
        #             print('at {:s}: skipping it'.format(cluster_path))
        #         # Go to next iteration
        #         continue
        #
        #     # Define current cluster HMM model
        #     hmm_path = os.path.join(cluster_path, 'HMM')
        #     # Check if HMM path exists
        #     if not os.path.isfile(hmm_path):
        #         # Remove cluster subfolder
        #         os.remove(cluster_path)
        #         # Verbose log
        #         if verbose:
        #             print('Could not find HMM model file', end=' ')
        #             print('for cluster {:s}'.format(cluster_name), end=' ')
        #             print('at {:s}: skipping it'.format(hmm_path))
        #         # Go to next iteration
        #         continue
        #
        #     # Open HMM file in read mode
        #     with open(hmm_path, 'r') as hmm_file:
        #         # Go through file line by line
        #         for line in hmm_file:
        #             # Retrieve length (None if current line is not HMM length)
        #             hmm_len = HMM.get_length(line)
        #             # Case retrieved HMM length is bigger than current one
        #             if (hmm_len is not None) and (hmm_len > hmmlib_len):
        #                 # Set current HMM model length as HMM library length
        #                 hmmlib_len = hmm_len
        #             # Append line to HMM library file
        #             hmmlib_file.write(line)
        #
        #     # Add current cluster name to the ones correctly done
        #     in_hmmlib.add(cluster_name)
        #
        # # Close HMM library file
        # hmmlib_file.close()

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Return HMM library path, number of input HMM models and maximum model  size
        return library_path, library_len, library_width

    # Search HMM against a target dataset
    def search_hmm_dataset(self, model_path, dataset_chunks, e_value, z_score, n_cores, client, verbose=False):
        """ Search HMM model/library against target dataset

        Searches a HMM model or a HMM library (which are handled in the same
        way by the search script) against a target dataset, hence against a
        list of chunks, using hmmsearch in a distributed manner.

        Returns two lists: one for TBLOUT files (i.e. sequence hits) and one for
        DOMTBLOUT files (i.e. domain hits). Notice that returned lists contain
        paths to tabular files which will remain into temporary directory if not
        manually deleted.

        Args
        hmm_path (list)         Path to searched HMM model/library
        dataset_chunks (list)   List of target dataset chunks
        e_value (float)         E-value to define a match significant
        z_score (int)           Z-score defining whole dataset length
        n_cores (int)           Number of cores available per job
        client (Client)         Client used to spawn jobs on cluster
        verbose (bool)          Whether to print verbose log or not

        Return
        (str)                   Path to directory holdinh search results

        Raise
        (CalledProcessError)    In case hmmsearch could not run
        """
        # Run futures safely
        try:

            # Intialize timers
            time_beg, time_end = time(), 0.0
            # Verbose
            if verbose:
                print('Searching HMM model', end=' ')
                print('at {:s}'.format(model_path), end=' ')
                print('against {:d}'.format(len(dataset_chunks)), end=' ')
                print('dataset chunks...', end=' ')

            # Set root path as the same folder holding model
            root_path = os.path.dirname(model_path)
            # Define temporary folder holding DOMTBLOUT and TBLOUT results
            results_path = mkdtemp(prefix='hmmsearch', dir=root_path)

            # Initialize futures
            futures = dict()
            # Loop through each target dataset
            for chunk in dataset_chunks:
                # Get chunk path
                chunk_path = chunk.path
                # Define current output file name
                out_name, _ = os.path.splitext(os.path.basename(chunk_path))
                # Define current output file path
                out_path = os.path.join(results_path, out_name)
                # Submit search function
                futures[out_name] = client.submit(
                    func=self.search_hmm_chunk,
                    model_path=model_path,
                    chunk_path=chunk_path,
                    out_path=out_path,
                    e_value=e_value,
                    z_score=z_score,
                    n_cores=n_cores,
                    hmmsearch=self.hmmsearch
                )
            # Retrieve results (tuples containing paths to tabluar output files)
            client.gather(futures)

            # Verbose
            if verbose:
                # Update timer
                time_end = time()
                # Show execution time
                print('done in {:.0f} seconds'.format(time_end - time_beg))

            # Return results directory
            return results_path

        except CalledProcessError as err:
            # Verbose
            if verbose:
                # Print error STDERR
                print(' '.join([
                    'Could not run distributed hmmsearch:',
                    'returned code {}:\n'.format(err.returncode),
                    err.stdout.strip()
                ]))
            # Raise error interrupting pipeline
            raise

    @staticmethod
    def search_hmm_chunk(model_path, chunk_path, out_path, e_value, z_score, n_cores, hmmsearch):
        """ Search HMM against a single dataset (chunk)

        Args
        model_path (str)        Path to input HMM model
        chunk_path (str)        Path to dataset chunk serach HMM against
        out_path (str)          Path to search results, suffix like `.tblout`
                                ore `.domtblout` will be added
        e_value (float)         E-value to define a match significant
        z_score (int)           Z-score defining whole dataset length
        n_cores (int)           Number of cores allocalble to run hmmsearch
        hmmsearch (HMMSearch)   Instance of hmmsearch object
        """
        # Define if input chunk is gzipped
        chunk_gzip = is_gzip(chunk_path)
        # Case input file is compressed
        if chunk_gzip:
            # Define chunk name (remove directories)
            chunk_name = os.path.basename(chunk_path)
            # Remove .gz extension
            chunk_name = re.sub(r'\.gz$', '', chunk_name)
            # Separate file extension from chunk name
            chunk_name, chunk_ext = os.path.splitext(chunk_name)

            # Define output directory path
            out_dir = os.path.dirname(out_path)
            # Make temporary file in output directory
            _, unzip_path = mkstemp(dir=out_dir, suffix=chunk_ext)

            # Unzip target path and overwrite given target dataset
            chunk_path = gunzip(in_path=chunk_path, out_path=unzip_path)

        # Initialize tblout file path
        tblout_path = out_path + '.tblout'
        # Initialize domtblout file path
        domtblout_path = out_path + '.domtblout'

        # Safely run hmmsearch
        try:

            # Run hmmsearch against given chunk
            hmmsearch.run(
                # Path to searched HMM model/library file
                hmm_path=model_path,
                # Path to sequences target dataset
                target_path=chunk_path,
                # # Path to search output
                # out_path=out_path,
                # Path to per-sequence output table
                tblout_path=tblout_path,
                # Path to per-domain output table
                domtblout_path=domtblout_path,
                # # Path to per-pfam-domain output table
                # pfamtblout_path=pfamtblout_path,
                # Define number of CPUs allocated
                num_cpus=n_cores,
                # Set e-value
                seq_e=e_value,
                # Set z-score
                seq_z=z_score
            )

            # Remove unzipped temporary file
            if chunk_gzip:
                os.remove(chunk_path)

        # Catch exception
        except CalledProcessError:
            # Remove unzipped temporary file
            if chunk_gzip:
                os.remove(chunk_path)
            # Raise the error interrupting pipeline
            raise

    @staticmethod
    def parse_search_scores(tblout_path, tblout_type='sequences', e_value=0.01, min_bits=25.0, max_bits=999999.99):
        # Initialize TBLOUT parser class
        tblout_cls = None
        # Case TBLOUT for sequences hits
        if tblout_type == 'sequences':
            tblout_cls = Tblout
        # Case DOMTBLOUT for domains hits
        elif tblout_type == 'domains':
            tblout_cls = Domtblout
        # Otherwise, raise error
        else:
            raise ValueError(' '.join([
                'Error: output table type can be "domains" or "sequences",',
                '{} set instead'.format(tblout_type)
            ]))

        # Read table
        tblout = tblout_cls(tblout_path)
        # Return hits
        return tblout.get_scores(e_value=e_value, min_bits=min_bits, max_bits=max_bits)

    @staticmethod
    def parse_search_hits(tblout_path, tblout_type='sequences', scores=dict()):
        # Initialize TBLOUT parser class
        tblout_cls = None
        # Case TBLOUT for sequences hits
        if tblout_type == 'sequences':
            tblout_cls = Tblout
        # Case DOMTBLOUT for domains hits
        elif tblout_type == 'domains':
            tblout_cls = Domtblout
        # Otherwise, raise error
        else:
            raise ValueError(' '.join([
                'Error: output table type can be "domains" or "sequences",',
                '{} set instead'.format(tblout_type)
            ]))

        # Read table
        tblout = tblout_cls(tblout_path)
        # Return hits
        return tblout.get_hits(scores)

    # Retrieve HMM hits from a list of tabular output files
    def parse_hmm_hits(self, tblout_paths, domtblout_paths, e_value, client):
        """ Retrieve sequence hits from list tabular output files

        Args
        tblout_paths (list)         List of strings representing paths to
                                    tabular TBLOUT files
        domtblout_paths (list)      List of strings representing paths to
                                    tabular DOMTBLOUT files
        e_value (float)             E-value threshold for determining scores

        Return
        (dict)      Dictionary mapping cluster name to sequences bit-score thresholds
        (dict)      Dictionary mapping cluster name to domains bit-score thresholds
        (dict)      Dictionary mapping cluster name to sequences hits
        (dict)      Dictionary mapping cluster name to domains hits
        """
        # Initialize sequence scores
        sequence_scores = dict()
        # Loop through each sequence TBLOUT file
        for tblout_path in tblout_paths:
            # Compute scores on current chunk
            scores = self.parse_search_scores(
                tblout_type='sequences',
                tblout_path=tblout_path,
                e_value=e_value,
                min_bits=25.0,
                max_bits=999999.99
            )
            # Update sequence scores
            sequence_scores = merge_scores(sequence_scores, scores)

        # Initialize sequence hits
        sequence_hits = dict()
        # Loop through each TBLOUT file
        for tblout_path in tblout_paths:
            # Retrieve sequence hits
            hits = self.parse_search_hits(
                tblout_type='sequences',
                tblout_path=tblout_path,
                scores=sequence_scores
            )
            # Update sequence hits
            sequence_hits = merge_hits(sequence_hits, hits)

        # Initialize domain scores
        domain_scores = dict()
        # Loop through each sequence DOMTBLOUT file
        for domtblout_path in domtblout_paths:
            # Retrieve domain scores
            scores = self.parse_search_scores(
                tblout_type='domains',
                tblout_path=domtblout_path,
                e_value=e_value,
                min_bits=25.0,
                max_bits=999999.99
            )
            # Update domain scores
            domain_scores = merge_scores(domain_scores, scores)

        # Initialize domain hits dict
        domain_hits = dict()
        # Loop through each DOMTBLOUT file
        for domtblout_path in domtblout_paths:
            # Retrieve hits
            hits = self.parse_search_hits(
                tblout_type='domains',
                tblout_path=domtblout_path,
                scores=domain_scores
            )
            # Update domain hits
            domain_hits = merge_hits(domain_hits, hits)

        # Return sequence scores and sequence hits
        return sequence_scores, domain_scores, sequence_hits, domain_hits

    # Search multiple HMM against UniProt
    def search_hmm_uniprot(self, clusters_path, uniprot_path, uniprot_shape, min_jobs=1, max_jobs=100, verbose=False):
        """ Search multiple HMM models against UniProt

        First, take all HMM available from clusters directory: if one cluster
        does not have HMM, discard it, otherwise append HMM to a unique HMM
        library.

        Then, search HMM library against cluster chunks and retrieve
        significant sequence accession numbers. These accession numbers will
        then be used to define wether a cluster has intersection with UniProt
        or not.

        Clusters which have at least one intersection with UniProt will then be
        removed and put into uniprot/ directory.
        """
        # Make HMM library
        library_path, library_len, library_width = self.make_hmm_library(
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Get maximum number of cores
        max_cores = self.scheduler.max_cores
        # Get maximum available memory per job
        max_memory = as_bytes(self.scheduler.max_memory)

        # Verbose log
        if verbose:
            print('Maximum resources allowed per job are', end=' ')
            print('{:d} cores and'.format(max_cores), end=' ')
            print('{:d} GB of memory'.format(math.floor(max_memory / 1e09)))

        # Get minimum required memory per CPU
        req_memory = HMMSearch.get_memory(
            model_len=library_width,
            longest_seq=uniprot_shape[1]
        )
        # Define maximum number of CPUs (cores) allocable to run hmmsearch
        num_cores = HMMSearch.get_cpus(
            model_len=library_width,
            max_memory=max_memory,
            max_cpus=max_cores,
            longest_seq=uniprot_shape[1]
        )
        # Check if memory is sufficient
        if not num_cores:
            # Raise new memory exception
            raise MemoryError(' '.join([
                'Unable to start HMM search against UniProt dataset:',
                'minimum required memory is {:d} GB'.format(math.ceil(req_memory / 1e09)),
                'while maximum allocable is {:d} GB'.format(math.floor(max_memory / 1e09))
            ]))

        # Verbose
        if verbose:
            print('Making HMM search against UniProt dataset', end=' ')
            print('on {:d} cores with'.format(num_cores), end=' ')
            print('{:d} GB of memory each'.format(math.floor(max_memory / num_cores / 1e09)))

        # Refine max memory
        max_memory = self.scheduler.max_memory
        # Require maximum resources
        with self.scheduler.adapt(min_jobs, max_jobs, cores=num_cores, memory=max_memory):

            # Make search against UniProt (retrieve results directory)
            results_path = self.search_hmm_dataset(
                model_path=library_path,
                dataset_chunks=self.uniprot,
                e_value=1000,
                z_score=uniprot_shape[0],
                n_cores=num_cores,
                client=self.scheduler.client,
                verbose=verbose
            )

        # Require default resources
        with self.scheduler.adapt(min_jobs, max_jobs):
            # Define path to tblout files
            tblout_paths = glob(os.path.join(results_path, '*.tblout'))
            # Define path to domtblout files
            domtblout_paths = glob(os.path.join(results_path, '*.domtblout'))
            # Parse search hits, retrieve significat scores and hits
            sequence_scores, sequence_hits, domain_scores, domain_hits = self.parse_hmm_hits(
                tblout_paths=tblout_paths,
                domtblout_paths=domtblout_paths,
                e_value=self.uniprot_e_value,
                client=self.scheduler.client
            )

        # Remove library (now useless)
        os.unlink(library_path)
        # Remove results directory (now useless)
        shutil.rmtree(results_path)

        # Define names of clusters matching UniProt
        in_uniprot = set(sequence_hits.keys()) & set(domain_hits.keys())
        # Loop through each cluster having a match in UniProt
        for cluster_name in in_uniprot:
            # Define cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Case cluster path exists
            if os.path.isdir(cluster_path):
                # Discard cluster
                shutil.move(cluster_path, os.path.join(uniprot_path, cluster_name))

        # Verbose
        if verbose:
            # Show number of clusters discarded
            print('{:d}'.format(len(in_uniprot)), end=' ')
            print('clusters discarded due to UniProt hits')

    # Search multiple HMM against MGnifam
    def search_hmm_mgnifam(self, clusters_path, mgnifam_shape, author_name, min_jobs=1, max_jobs=100, verbose=False):
        """ Search multiple HMM models against MGnifam

        First, take all HMM available from clusters directory: if one cluster
        does not have HMM, discard it, otherwise append HMM to a unique HMM
        library.

        Then, search HMM library against cluster chunks and retrieve
        significant sequence accession numbers. These accession numbers will
        then be used to define cluster members and retrieve sequences out from
        MGnifam fasta dataset.

        Afterwards, alignments among cluster members and the HMM model itself
        are carried out through hmmalign script.

        Raise
        (MemoryError)       In case given memory does not allow to process
                            hmmsearch even with just one core per job
        """
        # Make HMM library
        library_path, library_len, library_width = self.make_hmm_library(
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Get maximum number of cores
        max_cores = self.scheduler.max_cores
        # Get maximum available memory per job
        max_memory = as_bytes(self.scheduler.max_memory)

        # Verbose log
        if verbose:
            print('Maximum resources allowed per job are', end=' ')
            print('{:d} cores and'.format(max_cores), end=' ')
            print('{:d} GB of memory'.format(math.floor(max_memory / 1e09)))

        # Get minimum required memory per CPU
        req_memory = HMMSearch.get_memory(
            model_len=library_width,
            longest_seq=mgnifam_shape[1]
        )
        # Define maximum number of CPUs (cores) allocable to run hmmsearch
        num_cores = HMMSearch.get_cpus(
            model_len=library_width,
            max_memory=max_memory,
            max_cpus=max_cores,
            longest_seq=mgnifam_shape[1]
        )
        # Check if memory is sufficient
        if not num_cores:
            # Raise new memory exception
            raise MemoryError(' '.join([
                'Unable to start HMM search against UniProt dataset:',
                'minimum required memory is {:d} GB'.format(math.ceil(req_memory / 1e09)),
                'while maximum allocable is {:d} GB'.format(math.floor(max_memory / 1e09))
            ]))

        # Verbose
        if verbose:
            print('Making HMM search against UniProt dataset', end=' ')
            print('on {:d} cores with'.format(num_cores), end=' ')
            print('{:d} GB of memory each'.format(math.floor(max_memory / num_cores / 1e09)))

        # Refine max memory
        max_memory = self.scheduler.max_memory
        # Require maximum resources
        with self.scheduler.adapt(min_jobs, max_jobs, cores=num_cores, memory=max_memory):
            # Make search against UniProt (retrieve results directory)
            results_path = self.search_hmm_dataset(
                model_path=library_path,
                dataset_chunks=self.mgnifam,
                e_value=1000,
                z_score=mgnifam_shape[0],
                n_cores=num_cores,
                client=self.scheduler.client,
                verbose=verbose
            )

        # Require default resources
        with self.scheduler.adapt(min_jobs, max_jobs):
            # Define path to tblout files
            tblout_paths = glob(os.path.join(results_path, '*.tblout'))
            # Define path to domtblout files
            domtblout_paths = glob(os.path.join(results_path, '*.domtblout'))
            # Parse search hits, retrieve significat scores and hits
            sequence_scores, sequence_hits, domain_scores, domain_hits = self.parse_hmm_hits(
                tblout_paths=tblout_paths,
                domtblout_paths=domtblout_paths,
                e_value=self.mgnifam_e_value,
                client=self.scheduler.client
            )

            # Initailize cluster members dictionary
            cluster_members = dict()
            # Initialize sequences accession numbers
            sequences_acc = set()
            # Loop through each cluster name (query name) in domain hits
            for cluster_name in domain_hits:
                # Initialize cluster entry
                cluster_members.setdefault(cluster_name, set())
                # Go through each hit for current cluster name
                for hit in domain_hits[cluster_name]:
                    # Retrieve sequence accession (target name)
                    acc = hit.get('target_name')
                    # Update sequences accession
                    sequences_acc.add(acc)
                    # Update cluster members
                    cluster_members.get(cluster_name).add(acc)

            # Retrieve sequences related to retrieved accession numbers
            fasta_sequences, _ = self.search_fasta_sequences(
                sequences_acc=sequences_acc,
                fasta_dataset=self.mgnifam,
                client=self.scheduler.client,
                verbose=verbose
            )

            # Align cluster sequences to HMM using hmmalign
            self.make_hmm_alignments(
                clusters_path=clusters_path,
                fasta_sequences=fasta_sequences,
                cluster_members=cluster_members,
                client=self.scheduler.client
            )

            # Make clusters DESC files
            self.make_clusters_desc(
                clusters_path=clusters_path,
                sequence_scores=sequence_scores,
                domain_scores=domain_scores,
                sequence_hits=sequence_hits,
                domain_hits=domain_hits,
                author_name=author_name,
                verbose=verbose
            )

    # Make HMM alignment
    def make_hmm_alignments(self, clusters_path, fasta_sequences, cluster_members, client, verbose=False):
        """ Align all HMM models to sequneces

        Args
        clusters_path (str)     Path where cluster folders can be found
        fasta_sequences (dict)  Dictionary mapping accession number to fasta
                                sequences
        cluster_members (dict)  Dictionary associating cluster names to
                                seuqnece accession numbers
        client(dask.Client)     Dask client used to submit distributed jobs
        verbose (bool)          Whether to show verbose log or not

        Raise
        (OSError)               In case files can not be created
        (CalledProcessError)    In case hmmalign process can not be run
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making HMM model alignments for', end=' ')
            print('{:d} clusters'.format(len(cluster_members)), end=' ')
            print('at {:s}...'.format(clusters_path), end=' ')

        # Initialize futures
        futures = dict()
        # Loop through each cluster path
        for cluster_path in glob(os.path.join(clusters_path, 'MGYP*')):
            # Define cluster name
            cluster_name = os.path.basename(cluster_path)

            # Define path to fasta file
            fasta_path = os.path.join(cluster_path, 'ALIGN.fa')
            # Open file in write mode
            with open(fasta_path, 'w') as fasta_file:
                # Loop through every sequence in current cluster
                for sequence_acc in cluster_members.get(cluster_name, set()):
                    # Retrieve each fasta sequence
                    fasta_entry = fasta_sequences.get(sequence_acc)
                    # Write fasta sequence
                    fasta_file.write(fasta_entry + '\n')

            # Make alignment
            futures[cluster_name] = client.submit(
                func=self.align_hmm_cluster,
                cluster_path=cluster_path,
                hmmalign=self.hmmalign
            )
        # Gather all futures
        client.gather(futures)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    @staticmethod
    def align_hmm_cluster(cluster_path, hmmalign):
        """ Align HMM to cluster sequneces

        Raise
        (FileNotFoundError)     In case one of the input files cannot be opened
        (CalledProcessError)    In case hmmalign did not run properly
        """
        # Define cluster name
        cluster_name = os.path.basename(cluster_path)

        # Define input HMM model path
        model_path = os.path.join(cluster_path, 'HMM.model')
        # Define path to input fasta file
        fasta_path = os.path.join(cluster_path, 'ALIGN.fa')
        # Define output alignment path (.sto stockholm)
        sto_path = os.path.join(cluster_path, 'ALIGN.sto')
        # Define output  alignment path (.aln)
        aln_path = os.path.join(cluster_path, 'ALIGN.aln')

        # Check if HMM model file path is valid
        if not os.path.isfile(model_path):
            # Raise exception
            raise FileNotFoundError(' '.join([
                'Unable to run hmmalign script:',
                'HMM model for cluster {:s}'.format(cluster_name),
                'not found at {:s}'.format(model_path)
            ]))

        # Check if input FASTA file path is valid
        if not os.path.isfile(fasta_path):
            # Raise exception
            raise FileNotFoundError(' '.join([
                'Unable to run hmmalign script:',
                'input FASTA file for cluster {:s}'.format(cluster_name),
                'not found at {:s}'.format(fasta_path)
            ]))

        # Run hmmalign script
        hmmalign.run(
            # Set path to input HMM model
            hmm_path=model_path,
            # Set path to input FASTA file
            fasta_path=fasta_path,
            # Set path to output alignment file
            out_path=sto_path
        )
        # Load retrieved alignment file
        align = MSA.from_sto(sto_path)
        # Parse loaded alignment file
        align.to_aln(aln_path)

    # Make DESC files
    def make_clusters_desc(self, clusters_path, sequence_scores, domain_scores, sequence_hits, domain_hits, author_name='Anonymous', verbose=False):
        """ Make DESC files for each cluster

        Args
        clusters_path (str)         Root directory holding clusters folders
        sequence_scores (dict)      Dictionary mapping cluster names to
                                    bit scores thresholds for sequences
        domain_scores (dict)        Dictionary mapping cluster names to
                                    bit scores thresholds for domains
        sequence_hits (dict)        Dictionary mapping cluster names to
                                    list of significant sequence hits
        domain_hits (dict)          Dictionary mapping cluster names to
                                    list of significant domain hits
        author_name (str)           Name of the author responsible of running
                                    current build pipeline
        verbose (bool)              Whether to print out verbose log

        Raise
        (OSError)                   In case DESC file can not be made
        """
        # Define clusters paths iterator
        clusters_iter = glob(os.clusters_path, 'MGYP*')
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making DESC files for', end=' ')
            print('{:d} clusters...'.format(len(clusters_iter)), end=' ')

        # Loop through each cluster path
        for cluster_path in clusters_iter:
            # Get cluster name
            cluster_name = os.path.basename(cluster_path)
            # Make new mgnifam cluster DESC file
            cluster = Cluster(
                name=cluster_name,
                path=cluster_path,
                desc='Protein of unknown function (MGDUF)',
                auth=author_name,
                seq_scores=sequence_scores.get(cluster_name),
                dom_scores=domain_scores.get(cluster_name),
                type='Family'
            )

            # Define domain hits file path
            hits_path = os.path.join(cluster_path, 'HITS.tsv')
            # Write hits to .tsv file
            to_tsv(hits=domain_hits, path=hits_path, sep=r'\t')

            # Store to file
            cluster.to_desc()

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_beg - time_end))


# Utility: retrieve cluster names from TSV files
def get_cluster_names(paths, shuffle=False):
    """ Retrieve list of cluster names

    Given a list of file paths, opens them iteratively and extracts cluster
    names. Optionally, shuffles results to avoid undesired correlations.

    Args
    paths (list)            List of input file paths
    shuffle (bool)          Wether to shuffle cluster names before returning
                            them or not

    Return
    (list)                  List of cluster names as strings

    Raise
    (FileNotFoundError)     In case at least one file provided as input is not
                            found at given path
    """
    # Initialize outout cluster names list
    cluster_names = list()
    # Loop through each input path
    for i in range(len(paths)):
        # Open file
        file = open(paths[i], 'r')
        # Iterate through each file line
        for line in file:
            # Check if current line matches expected format
            match = re.search(r'^(\S+)', line)
            # Case current line does not match expected format
            if not match:
                # Skip line, go to next iteration
                continue
            # Otherwise, store cluster name
            cluster_names.append(match.group(1))
        # Close file
        file.close()
    # Check if shuffle is active
    if shuffle:
        # Shuffle output list before returning it
        random.shuffle(cluster_names)
    # Return cluster names list
    return cluster_names


# Utility: get environmental variables
def parse_env_json(path):
    # Case environment .json file path does not exist
    if not os.path.isfile(path):
        # Raise error
        raise FileNotFoundError(' '.join([
            'could not open environmental variables file',
            'at %s: file not found'.format(path)
        ]))

    # Initialize environmental variables dictionary
    env = dict()
    # Open environmental variables file
    with open(path, 'r') as file:
        # Load variables dictionary from file
        env = json.load(file)
        # Loop through every key in retrieved dictionary
        for key, value in env.items():
            # Case calue is not a list
            if isinstance(value, list):
                # Join list as string, using double dots as separators
                env[key] = ':'.join(value)

    # Return environmental variables dictionary
    return env


# Run batch pipeline
if __name__ == '__main__':

    # Define project root path
    ROOT_PATH = os.path.dirname(__file__) + '/../..'
    # Define path to data folder
    DATA_PATH = ROOT_PATH + '/data'
    # Define path to LinClust file(s)
    LINCLUST_PATH = DATA_PATH + '/clusters/chunk*.tsv.gz'
    # Define path to MGnifam file(s)
    MGNIFAM_PATH = DATA_PATH + '/mgnify/chunk*.fa.gz'
    # Define path to UniProt file(s)
    UNIPROT_PATH = DATA_PATH + '/uniprot/chunk*.fa.gz'

    # Argparse
    parser = argparse.ArgumentParser(description='Build MGnifam clusters')
    # Define input path(s)
    parser.add_argument(
        '-i', '--in_path', nargs='+', type=str, required=True,
        help='Path to input .tsv file(s), cluster name is first column'
    )
    # Define output path
    parser.add_argument(
        '-o', '--out_path', type=str, required=True,
        help='Path to output directory'
    )
    # Define author name
    parser.add_argument(
        '-a', '--author_name', type=str, default='Anonymous',
        help='Name of the user running MGnifam build pipeline'
    )
    # Wether to shuffle input or not
    parser.add_argument(
        '--shuffle', type=int, default=0,
        help='Whether to shuffle input cluster names or not'
    )
    # Define batch size
    parser.add_argument(
        '--batch_size', type=int, default=100,
        help='Number of clusters to make at each iteration'
    )
    # Define maximum number of clusters to make
    parser.add_argument(
        '--max_clusters', type=int, default=0,
        help='Maximum number of clusters to make'
    )
    # Define path to LinClust clusters file(s)
    parser.add_argument(
        '--linclust_path', nargs='+', type=str, default=glob(LINCLUST_PATH),
        help='Path to LinClust clusters file(s)'
    )
    # Define path to MGnifam file(s)
    parser.add_argument(
        '--mgnifam_path', nargs='+', type=str, default=glob(MGNIFAM_PATH),
        help='Path to MGnifam file(s)'
    )
    # Define MGnifam width (maximum sequence length)
    parser.add_argument(
        '--mgnifam_width', type=int, required=False,
        help='Maximum sequence length in MGnifam dataset'
    )
    # Define MGnifam height (total number of sequences)
    parser.add_argument(
        '--mgnifam_height', type=int, required=False,
        help='Total number of sequences in MGnifam dataset'
    )
    # Define path to UniProt file(s)
    parser.add_argument(
        '--uniprot_path', nargs='+', type=str, default=glob(UNIPROT_PATH),
        help='Path to UniProt file(s)'
    )
    # Define UniProt width (maximum sequence length)
    parser.add_argument(
        '--uniprot_width', type=int, required=False,
        help='Maximum sequence length in UniProt dataset'
    )
    # Define UniProt height (total number of sequences)
    parser.add_argument(
        '--uniprot_height', type=int, required=False,
        help='Total number of sequences in UniProt dataset'
    )
    # Define MobiDB executable
    parser.add_argument(
        '--mobidb_cmd', nargs='+', type=str, default=['mobidb_lite.py'],
        help='MobiDB Lite disorder predictor executable'
    )
    # Define Muscle executable
    parser.add_argument(
        '--muscle_cmd', nargs='+', type=str, default=['muscle'],
        help='Muscle multiple sequence aligner executable'
    )
    # Define hmmsearch executable
    parser.add_argument(
        '--hmmsearch_cmd', nargs='+', type=str, default=['hmmsearch'],
        help='HMMER3 search executable'
    )
    # Define hmmbuild executable
    parser.add_argument(
        '--hmmbuild_cmd', nargs='+', type=str, default=['hmmbuild'],
        help='HMMER3 build executable'
    )
    # Define hmmalign executable
    parser.add_argument(
        '--hmmalign_cmd', nargs='+', type=str, default=['hmmalign'],
        help='HMMER3 align executable'
    )
    # Whether to print verbose output
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Print verbose output'
    )
    # Define E-value significance threshold for both UniProt and MGnifam
    parser.add_argument(
        '-e', '--e_value', type=float, default=0.01,
        help='E-value threhsold for both UniProt and MGnifam comparisons'
    )
    # Define environmental variables .json file
    parser.add_argument(
        '--env_path', type=str, required=False,
        help='Path to .json file holding environmental variables'
    )
    # Define scheduler options
    group = parser.add_argument_group('Scheduler options')
    # Define schduler type
    group.add_argument(
        '-s', '--scheduler_type', type=str, default='LSF',
        help='Type of scheduler to use to distribute parallel processes'
    )
    # Define minimum number of jobs
    group.add_argument(
        '-j', '--min_jobs', type=int, default=0,
        help='Minimum number of parallel processes to keep alive'
    )
    # Define maximum number of jobs
    group.add_argument(
        '-J', '--max_jobs', type=int, default=100,
        help='Maximum number of parallel processes to keep alive'
    )
    # Define minimum number of cores
    group.add_argument(
        '-c', '--min_cores', type=int, default=1,
        help='Minimum number of cores to use per process'
    )
    # Define maximum number of cores
    group.add_argument(
        '-C', '--max_cores', type=int, default=1,
        help='Maximum number of cores to use per process'
    )
    # Define minimum memory allocable per job
    group.add_argument(
        '-m', '--min_memory', type=str, default='1 GB',
        help='Minimum memory allocable per process'
    )
    # Define maximum memory allocable per job
    group.add_argument(
        '-M', '--max_memory', type=str, default='1 GB',
        help='Maximum memory allocable per process'
    )
    # Define walltime
    group.add_argument(
        '-W', '--walltime', type=str, default='02:00',
        help='How long can a process be kept alive'
    )
    # Retrieve arguments
    args = parser.parse_args()

    # Initialize scheduler
    scheduler = None
    # Case scheduler type is LSF
    if args.scheduler_type == 'LSF':
        # Define LSF scheduler
        scheduler = LSFScheduler(
            # Define cores boundaries
            min_cores=args.min_cores,
            max_cores=args.max_cores,
            # Define memory boundaries
            min_memory=args.min_memory,
            max_memory=args.max_memory,
            # Debug
            silence_logs='debug',
            # Define walltime
            walltime=args.walltime,
            # Define processes per job
            processes=1
        )
    # Case scheduler type is Local
    if args.scheduler_type == 'Local':
        # Define Local scheduler
        scheduler = LocalScheduler(
            # Define cores boundaries
            min_cores=args.min_cores,
            max_cores=args.max_cores,
            # Define memory boundaries
            min_memory=args.min_memory,
            max_memory=args.max_memory,
            # Define processes per job
            threads_per_worker=1
        )
    # Case no scheduler has been set
    if scheduler is None:
        # Raise new error
        raise ValueError(' '.join([
            'scheduler type can be one among `Local` or `LSF`',
            '%s has been chosen instead' % args.scheduler_type
        ]))

    # Get cluster names list
    cluster_names = get_cluster_names(
        paths=args.in_path,
        shuffle=args.shuffle
    )

    # Load environmental variables
    env = os.environ.copy()
    # Case environmental variable file is set
    if args.env_path:
        # Update environmental variables using given file file
        env = {**env, **parse_env_json(args.env_path)}

    # Try making release directory
    try:
        # Case output directory does not exist
        if not os.path.isdir(args.out_path):
            # Make output directory
            os.mkdir(args.out_path)

    except OSError:
        # Raise new error
        raise OSError(' '.join([
            'could not create output directory',
            'at %s:' % args.out_path,
            'permission denied'
        ]))

    # # Load dataset chunks
    # linclust_chunks = LinClust.from_list(args.linclust_path)
    # mgnifam_chunks = Fasta.from_list(args.mgnifam_path)
    # uniprot_chunks = Fasta.from_list(args.uniprot_path)

    # Define new build pipeline
    build = Build(
        # Set scheduler
        scheduler=scheduler,
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
        hmmbuild_cmd=args.hmmbuild_cmd,  # Path to hmmbuild script
        hmmsearch_cmd=args.hmmsearch_cmd,  # Path to hmmsearch script
        hmmalign_cmd=args.hmmalign_cmd,  # Path to hmmalign script
        # Environmental variables
        env=env
    )

    # Define number of clusters
    num_clusters = len(cluster_names)
    # Threshold number of clusters
    num_clusters = min(num_clusters, args.max_clusters)
    # Define batchs size
    batch_size = args.batch_size
    # Loop through batches of input cluster names
    for i in range(0, num_clusters, batch_size):
        # Define last index of current batch (excluded)
        j = min(num_clusters, i + batch_size)
        # Define batch index
        batch_index = str(i // batch_size)
        # Define directory for current batch
        batch_path = os.path.join(args.out_path, 'batch', batch_index)

        # Try making output directory
        try:
            # Make directory for current batch
            os.makedirs(batch_path, exist_ok=True)

        # Intercept exception
        except OSError:
            # Raise exception
            raise OSError(' '.join([
                'could not make batch directory',
                'at %s:'.format(batch_path),
                'permission denied'
            ]))

        # Define batch of cluster names
        batch_names = cluster_names[i:j]
        # Run pipeline for current batch of cluster names
        build(
            author_name=args.author_name,
            cluster_names=batch_names,
            clusters_path=batch_path,
            min_jobs=args.min_jobs,
            max_jobs=args.max_jobs,
            verbose=bool(args.verbose)
        )
