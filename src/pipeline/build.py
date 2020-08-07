# Dependencies
from src.pipeline.pipeline import Pipeline
from src.msa.transform import Compose, OccupancyTrim, OccupancyFilter
from src.msa.msa import MSA, Muscle
from src.hmm.hmmer import HMMBuild, HMMAlign
from src.hmm.hmmsearch import HMMSearch, Tblout, Domtblout
from src.hmm.hmmsearch import merge_scores, merge_hits
from src.hmm.hmm import HMM
from src.dataset import LinClust
from src.dataset import Fasta
from src.disorder import MobiDbLite
from src.utils import benchmark, is_gzip, gunzip, as_bytes
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
from glob import iglob, glob
from time import time
import shutil
import random
import math
import json
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
        self, cluster_type, cluster_kwargs,
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
        hmm_build_cmd=['hmmbuild'],  # Path to hmmbuild script
        hmm_align_cmd=['hmmalign'],  # Path to hmmalign script
        hmm_search_cmd=['hmmsearch'],  # Path to hmmsearch script
        # Environmental variables
        env=os.environ.copy()
    ):
        # Call parent constructor
        super().__init__(cluster_type, cluster_kwargs)
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
        self.hmm_build = HMMBuild(cmd=hmm_build_cmd, env=env)
        self.hmm_align = HMMAlign(cmd=hmm_align_cmd, env=env)
        self.hmm_search = HMMSearch(cmd=hmm_search_cmd, env=env)
        # Store environmental variables
        self.env = env

    # Run the pipeline
    def run(self, cluster_names, clusters_path, min_jobs=1, max_jobs=100, verbose=False):
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

        # Initialize timers
        time_beg, time_end = time(), 0.0

        # Verbose log
        if verbose:
            print('Running BATCH pipeline', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('at {:s}'.format(clusters_path))

        # Initialize sub-directories
        bias_path, noaln_path, uniprot_path = self.init_subdir(clusters_path)

        # Initialize new cluster
        cluster, client = self.get_client(dict(cores=1))
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Fill cluster members dictionary
        time_run, cluster_members = benchmark(
            fn=self.search_cluster_members,
            cluster_names=cluster_names,
            client=client,
            verbose=False
        )

        # Check cluster members
        self.check_cluster_members(cluster_names, cluster_members)

        # Initialize set of retrieved sequence accession numbers
        sequences_acc = set()
        # Loop through each cluster
        for members_acc in cluster_members.values():
            # Update sequence accession numbers
            sequences_acc |= set(members_acc)

        # Get number of clusters
        num_clusters = len(cluster_members)
        # Get number of sequences
        num_sequences = sum([len(c) for c in cluster_members.values()])

        # Verbose log
        if verbose:
            print('Retrieved {:d} accession numbers'.format(num_sequences), end=' ')
            print('for {:d} clusters'.format(num_clusters), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Fill mgnifam dictionary
        time_run, (fasta_sequences, _) = benchmark(
            fn=self.search_fasta_sequences,
            sequences_acc=sequences_acc,
            fasta_dataset=self.mgnifam,
            client=client,
            verbose=verbose
        )

        # Check fasta sequences
        self.check_fasta_sequences(sequences_acc, fasta_sequences)

        # Define number of sequences
        num_sequences = len(fasta_sequences)

        # Verbose log
        if verbose:
            print('Retrieved {:d} FASTA'.format(num_sequences), end=' ')
            print('sequences in {:.0f} seconds'.format(time_run))

        # Check compositional bias
        time_run, (cluster_names, _) = benchmark(
            fn=self.check_comp_bias,
            clusters_path=clusters_path,
            bias_path=bias_path,
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            client=client,
            verbose=verbose
        )

        # Keep only clusters whose bias is below threshold
        cluster_members = {
            cluster_name: sequences_acc
            for cluster_name, sequences_acc
            in cluster_members.items()
            if cluster_name in set(cluster_names)
        }

        # Verbose log
        if verbose:
            print('Compositional bias computed', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Make SEED alignments
        time_run, _ = benchmark(
            fn=self.make_seed_alignments,
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            clusters_path=clusters_path,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print('SEED sequence alignments generated', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Trim SEED alignments
        time_run, cluster_names = benchmark(
            fn=self.trim_seed_alignments,
            cluster_names=set(cluster_members.keys()),
            clusters_path=clusters_path,
            noaln_path=noaln_path,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print('SEED sequence alignments trimmed', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Make HMM models
        time_run, cluster_names = benchmark(
            fn=self.make_hmm_models,
            cluster_names=cluster_names,
            clusters_path=clusters_path,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print('HMM models done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Check if UniProt shape hasn't been already set
        if not (self.uniprot_height and self.uniprot_width):
            # Retrieve UniProt size and it longest sequence length
            time_run, (_, self.uniprot_width, self.uniprot_height) = benchmark(
                fn=self.get_longest,
                target_dataset=self.uniprot,
                client=client
            )
            # Verbose log
            if verbose:
                print('UniProt shape retrieved in {:.0f} seconds'.format(time_run))
        # Define UniProt shape tuple
        uniprot_shape = (self.uniprot_height, self.uniprot_width)

        # Verbose log
        if verbose:
            print('UniProt dataset has shape', end=' ')
            print('{:d} x {:d}'.format(*uniprot_shape))

        # Check if MGnifam shape hasn't been already set
        if not (self.mgnifam_height and self.mgnifam_width):
            # Retrieve MGnifam size and it longest sequence length
            time_run, (_, self.mgnifam_width, self.mgnifam_height) = benchmark(
                fn=self.get_longest,
                target_dataset=self.mgnifam,
                client=client
            )
            # Verbose log
            if verbose:
                print('MGnifam shape retrieved in {:.0f} seconds'.format(time_run))
        # Store uniprot shape
        mgnifam_shape = (self.mgnifam_height, self.mgnifam_width)

        # Verbose log
        if verbose:
            print('MGnifam dataset has shape', end=' ')
            print('{:d} x {:d}'.format(*mgnifam_shape))

        # Close previous client
        self.close_client(cluster, client)

        # Run search against UniProt sub-pipeline
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
            min_jobs=min_jobs,
            max_jobs=max_jobs,
            verbose=verbose
        )

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Define clusters path iterator
        cluster_iter = glob(os.path.join(clusters_path, 'MGYP*'))
        # Define cluster names
        cluster_names = set(map(os.path.basename, cluster_iter))

        # Verbose log
        if verbose:
            print('Build pipeline done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.0f} seconds'.format(time_tot))

        # Return time taken to run and set of cluster names
        return time_tot, cluster_names

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
        sub_paths = ['bias', 'noaln', 'uniprot']
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
    def get_longest(self, target_dataset, client):
        """ Retrieve longest sequence and dataset length

        Args
        target_dataset (list)       List of target dataset chunks to go through
        client (Client)             Dask client used to distribute search

        Return
        (str)                       Longest sequence in FASTA dataset
        (int)                       Length of longest sequence in FASTA dataset
        (int)                       Dataset length
        """
        # Initailize list of futures
        futures = list()
        # Loop through each chunk in dataset
        for chunk in target_dataset:
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
        # Verbose log
        if verbose:
            print('Searching member sequences', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)))

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

        # Verbose log
        if verbose:
            print('Retrieved sequences', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)))

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
        # Verbose log
        if verbose:
            print('Searching {:d}'.format(len(sequences_acc)), end=' ')
            print('FASTA sequences')

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

        # Verbose log
        if verbose:
            print('Retrieved {:d} fasta'.format(len(fasta_sequences)), end=' ')
            print('entries from dataset of size {}'.format(fasta_length))

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

    # Check compositional bias
    def check_comp_bias(self, clusters_path, bias_path, cluster_members, fasta_sequences, client, verbose=False):
        """ Check compositional bias

        Given a cluster members dictionary and a fasta sequences dictionary,
        makes clusters dict(cluster members: fasta sequences). Then makes
        disorder predictions over the sequences: if rate of disordered residues
        is too high, cluster is discarded

        Args
        cluster_members (dict)      Dictionary mappung cluster names to cluster
                                    seqeunces accession numbers
        fasta_sequences (dict)      Dictionary mapping seqeunces accession
                                    numbers to fasta entries
        client (Client)             Dask client used to submit jobs
        verbose (bool)              Wether to print verbose output or not

        Return
        (dict)                      Input cluster entries above compositional
                                    bias threshold
        (dict)                      Input cluster entries below compositional
                                    bias threshold
        """
        # Verbose log
        if verbose:
            print('Computing compositional bias', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('and {:d} sequences'.format(len(fasta_sequences)))

        # Initialize futures dict(cluster name: compositional bias)
        futures = dict()
        # Loop through each cluster
        for cluster_name in cluster_members:

            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = dict()
            # Loop through each sequence accession in cluster
            for acc in cluster_members.get(cluster_name):
                # Retrieve fasta sequence
                cluster_sequences[acc] = fasta_sequences[acc]

            # Run distributed compositional bias computation
            futures[cluster_name] = client.submit(
                func=self.compute_comp_bias,
                fasta_sequences=cluster_sequences,
                pred_threshold=1,  # Set threshold for disorder prediction
                pred_inclusive=True,  # Set thresold inclusive
                mobidb=self.mobidb
            )

        # Retrieve results
        results = client.gather(futures)

        # Initialize set of clusters below and above threshold
        bias_below, bias_above = set(), set()
        # Retrieve compositional bias threshold
        threshold = self.comp_bias_threshold
        # Retrieve compositional bias threshold is inclusive or not
        inclusive = self.comp_bias_inclusive
        # Loop through each cluster name in results
        for cluster_name in results.keys():

            # Define current cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Encure that cluster path exists
            if not os.path.isdir(cluster_path):
                # Make new directory
                os.mkdir(cluster_path)

            # Retrieve cluster bias
            cluster_bias = results.get(cluster_name, 1.0)

            # Write bias score to file
            with open(os.path.join(cluster_path, 'BIAS'), 'w') as file:
                file.write(str(cluster_bias))

            # Define target path if bias is too high
            target_path = os.path.join(bias_path, cluster_name)
            # Case compositional bias threshold is inclusive
            if (cluster_bias >= threshold) and inclusive:
                # Move cluster subfolder into BIAS folder
                shutil.move(cluster_path, target_path)
                # Add cluster name to the above threshold set
                bias_above.add(cluster_name)
                # Go to next iteration
                continue
            # Case compositional bias threshold is exclusive
            if (cluster_bias > threshold) and not inclusive:
                # Move cluster subfolder into BIAS folder
                shutil.move(cluster_path, target_path)
                # Add cluster name to the above threshold set
                bias_above.add(cluster_name)
                # Go to next iteration
                continue

            # By default, add cluster name to the below threshold set
            bias_below.add(cluster_name)

        # Verbose log
        if verbose:
            print('Compositional bias computed for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('among which {:d} were below'.format(len(bias_below)), end=' ')
            print('threshold of {:.02f},'.format(self.comp_bias_threshold), end=' ')
            print('while {:d} were above it'.format(len(bias_above)))

        # Return both below and above threshold sets
        return bias_below, bias_above

    @staticmethod
    def compute_comp_bias(fasta_sequences, mobidb, pred_threshold=1, pred_inclusive=True):
        """ Compute compositional bias for a cluster

        Given a sequences dict(seqeunce acc: fasta entry), run MobiDbLite
        predictor against them and produces a disordered prediction for each
        residue. A residue is considered to be disordere if it has at least
        <pred_threshold> disordered predictions.

        Once disordered predictiona are computed, compositional bias is
        calculated as rate of predicted disordered residues over number of
        total residues.

        Args
        seqeunces (dict)        Dictionary mapping sequence accession to actual
                                FASTA entry
        mobidb (MobiDbLite)     Instance of MobiDB Lite interface
        pred_threshold (int)    Number of times a residue must be predicted
                                disordered to be considered disordered.
        pred_inclusive (bool)   Wether the prediction threshold should be
                                inclusive or not

        Return
        (float)                 Compositional bias score
        """
        # Make predictions running MobiDB Lite
        predictions = mobidb.run(fasta_sequences=fasta_sequences)
        # Apply threshold over predictions
        predictions = mobidb.apply_threshold(
            sequences=predictions,
            threshold=pred_threshold,
            inclusive=pred_inclusive
        )
        # Compute and return compositional bias
        return mobidb.compute_bias(predictions)

    # Make SEED alignments
    def make_seed_alignments(self, cluster_members, fasta_sequences, clusters_path, client, verbose=False):
        """ Make SEED alignments

        Given a cluster memebers dictionary associating sequences accession
        numbers to a cluster name and a fasta sequences dictionary associating
        to each sequence accession number its fasta entry, SEED alignment is
        produced through Muscle algorithm and stored at given cluster directory.

        Args
        cluster_members (dict)      Dictionary associating cluster names (keys)
                                    to its member sequences accession numbers
                                    (values)
        fasta_seqeunces (dict)      Dictionary associating sequence accession
                                    numbers (keys) to fasta entries (values)
        clusters_path (str)         Path where cluster folders can be found
        client (str)                Client used to schedule jobs
        verbose (bool)              Wether to show verbose log or not
        """
        # Verbose log
        if verbose:
            print('Making SEED alignments', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('at {:s}'.format(clusters_path))

        # Initialize futures
        futures = list()
        # Loop through each cluster
        for cluster_name in cluster_members:

            # Define current cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Ensure that cluster directory exists
            if not os.path.exists(cluster_path):
                # Try to make one
                os.mkdir(cluster_path)

            # Define raw alignment path
            alignment_path = os.path.join(cluster_path, 'SEED_raw')

            # Initialize cluster sequences
            cluster_sequences = list()
            # Loop through each sequence accession in cluster
            for sequence_acc in cluster_members.get(cluster_name):
                # Save fasta sequence
                cluster_sequences.append(fasta_sequences[sequence_acc])

            # Make SEED alignment in cluster directory
            future = client.submit(
                func=self.make_muscle_alignment,
                fasta_sequences=cluster_sequences,
                alignment_path=alignment_path,
                muscle=self.muscle,
            )

            # Store future
            futures.append(future)

        # Retrieve results
        client.gather(futures)

    @staticmethod
    def make_muscle_alignment(fasta_sequences, alignment_path, muscle):
        """ Use MUSCLE to make a multiple sequence alignment

        Args
        fasta_sequences (list)      List of fasta sequences
        alignment_path (str)        Path where alignment will be moved once
                                    done
        muscle (Muscle)             Instance of Muscle interface object
        """
        # Make multiple sequence alignment
        alignment = muscle.run(
            sequences=fasta_sequences,
            acc_regex='>(\S+)',
            verbose=False
        )
        # Store alignment to disk
        alignment.to_aln(alignment_path)

    # Automated SEED alignments trimming
    def trim_seed_alignments(self, cluster_names, clusters_path, noaln_path, verbose=False):
        """ Trim SEED alignments of clusters in given directory

        Goes into each cluster (directory) in given clusters path, then tries to
        load alignment file and to trim it. If alignment cannot be loaded, then
        it is skipped and moved in noaln/ folder.

        Args
        cluster_names (set)     Set of cluster names whose SEED has to be
                                automatically trimmed
        clusters_path (str)     Path to directory holding cluster directories
        noaln_path (str)        Path where to put clusters whose trimmed SEED
                                alignment does not fulfill requirements
        verbose (bool)          Whether to print verbose log or not

        Return
        (set)                   Set of cluster names which fulfill the
                                trimmed SEED alignment requirements
        """

        # Verbose log
        if verbose:
            print('Trimming SEED alignments', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('at {:s}'.format(clusters_path))

        # Define set of cluster names whose SEED alignment has been trimmed
        done_clusters = set()

        # Loop through each cluster name
        for cluster_name in cluster_names:

            # Define current cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Define path to raw alignment
            raw_seed_path = os.path.join(cluster_path, 'SEED_raw')
            # Define path to new alignment
            new_seed_path = os.path.join(cluster_path, 'SEED')

            # Ensure that the curreny cluster directory exists
            if not os.path.exists(raw_seed_path):
                # Verbose log
                if verbose:
                    print('Could not find SEED alignment', end=' ')
                    print('for cluster {:s}'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(cluster_path))
                # Skip iteration
                continue

            # Load SEED alignment from file
            raw_seed_aln = MSA.from_aln(raw_seed_path)
            # Trim SEED alignment
            new_seed_aln = self.seed_transform(raw_seed_aln)

            # Retrieve new SEED alignment shape
            n, m = new_seed_aln.aln.shape
            # Case new SEED does not fulfill requirements
            if (n < self.seed_min_height) or (m < self.seed_min_width):

                # Verbose log
                if verbose:
                    print('Trimmed SEED alignment', end=' ')
                    print('for cluster {:s}'.format(cluster_name), end=' ')
                    print('has shape {:d} x {:d}:'.format(n, m), end=' ')
                    print('skipping it')

                # Move cluster to noaln path
                shutil.move(cluster_path, os.path.join(noaln_path, cluster_name))

                # Skip iteration
                continue

            # Save trimmed SEED alignment to file
            new_seed_aln.to_aln(new_seed_path)

            # Add cluster name to valid cluster names
            done_clusters.add(cluster_name)

        # Return done cluster names
        return done_clusters

    # Make HMM models
    def make_hmm_models(self, cluster_names, clusters_path, client, verbose=False):
        """ Build HMM models

        For each cluster in given clusters path, create a new HMM model through
        hmmbuild script, then store it in the associated cluster directory.

        Args
        cluster_names (set)     Set of cluster names whose HMM has to be done
        clusters_path (str)     Path where cluster directories can be found
        client (Client)         Client used to spawn jobs on cluster
        verbose (bool)          Whether to print verbose log or not

        Return
        (set)                   Set of clusters which ran successfully
        """
        # Verbose log
        if verbose:
            print('Making HMM models', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)))

        # Initialize futures
        futures = dict()
        # Loop through each given cluster name
        for cluster_name in cluster_names:

            # Define cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Define SEED alignment path
            seed_path = os.path.join(cluster_path, 'SEED')
            # Define HMM path
            hmm_path = os.path.join(cluster_path, 'HMM')

            # Ensure that cluster directory exists
            if not os.path.exists(cluster_path):
                # Verbose log
                if verbose:
                    print('Cluster {:s} not found'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(cluster_path))
                # Skip iteration
                continue

            # Ensure that cluster SEED alignment exists
            if not os.path.exists(seed_path):
                # Verbose log
                if verbose:
                    print('SEED alignment for cluster', end=' ')
                    print('{:s} not found'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(cluster_path))
                # Skip iteration
                continue

            # Submit new job making HMM model
            futures[cluster_name] = client.submit(
                func=self.hmm_build.run,
                msa_path=seed_path,
                hmm_path=hmm_path,
                name=cluster_name
            )

        # Retrieve results
        results = client.gather([*futures.values()])
        # Return set of cluster names which ran successfully
        return set(futures.keys())

    # Make HMM library containing all HMM models
    def make_hmm_library(self, cluster_names, clusters_path, verbose=False):
        """ Concatenate every HMM mode in one HMM library

        If a cluster subfolder is not found: skip it
        If a cluster subfolder is found but inner HMM model is not found,
        remove whole cluster subfolder.

        Args
        cluster_names (set)     Set of cluster names whose HMM has to be done
        clusters_path (str)     Path where cluster directories can be found
        verbose (bool)          Whether to print verbose log or not

        Return
        (str)                   Path to HMM library
        (int)                   Length of longest HMM model in library
        (set)                   Set of cluster names whose HMM is in library
        """
        # Verbose log
        if verbose:
            print('Making HMM library', end=' ')
            print('out of {:d} clusters'.format(len(cluster_names)))

        # Initialize output HMM library path
        hmmlib_path = os.path.join(clusters_path, 'HMMLIB')
        # Open file in write mode
        hmmlib_file = open(hmmlib_path, 'w')

        # Initialize names of clusters loaded into HMM library
        in_hmmlib = set()
        # Initialize max HMM model model length in MM library
        hmmlib_len = 0
        # Loop through each cluster name in given cluster names set
        for cluster_name in cluster_names:

            # Define cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Check if cluster path exists
            if not os.path.isdir(cluster_path):
                # Verbose log
                if verbose:
                    print('Could not find subfolder', end=' ')
                    print('for cluster {:s}'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(cluster_path))
                # Go to next iteration
                continue

            # Define current cluster HMM model
            hmm_path = os.path.join(cluster_path, 'HMM')
            # Check if HMM path exists
            if not os.path.isfile(hmm_path):
                # Remove cluster subfolder
                os.remove(cluster_path)
                # Verbose log
                if verbose:
                    print('Could not find HMM model file', end=' ')
                    print('for cluster {:s}'.format(cluster_name), end=' ')
                    print('at {:s}: skipping it'.format(hmm_path))
                # Go to next iteration
                continue

            # Open HMM file in read mode
            with open(hmm_path, 'r') as hmm_file:
                # Go through file line by line
                for line in hmm_file:
                    # Retrieve length (None if current line is not HMM length)
                    hmm_len = HMM.get_length(line)
                    # Case retrieved HMM length is bigger than current one
                    if (hmm_len is not None) and (hmm_len > hmmlib_len):
                        # Set current HMM model length as HMM library length
                        hmmlib_len = hmm_len
                    # Append line to HMM library file
                    hmmlib_file.write(line)

            # Add current cluster name to the ones correctly done
            in_hmmlib.add(cluster_name)

        # Close HMM library file
        hmmlib_file.close()

        # Return HMM library path and length
        return hmmlib_path, hmmlib_len, in_hmmlib

    @staticmethod
    def search_hmm_chunk(hmm_path, chunk_path, e_value, z_score, n_cores, hmm_search):
        """ Search HMM against a single dataset (chunk)

        Args
        hmm_path (str)          Path to input HMM model
        chunk_path (str)        Path to dataset chunk serach HMM against
        e_value (float)         E-value to define a match significant
        z_score (int)           Z-score defining whole dataset length
        n_cores (int)           Number of cores allocalble to run hmmsearch
        hmm_search (HMMSearch)  Instance of hmmsearch object

        Return
        (str)                   Path where temporary file holding DOMTBLOUT
                                results is stored
        """

        # Define if input chunk is gzipped
        chunk_gzip = is_gzip(chunk_path)
        # Case input file is compressed
        if chunk_gzip:
            # Get file extension from chunk path
            _, chunk_ext = os.path.splitext(chunk_path)
            # Remove .gz extension
            chunk_ext = re.sub(r'\.gz$', '', chunk_ext)
            # Unzip target path and overwrite given target dataset
            chunk_path = gunzip(in_path=chunk_path, out_suffix=chunk_ext)

        # Initialize tblout temporary file path
        tblout_path = NamedTemporaryFile(delete=False, dir=os.path.dirname(hmm_path), suffix='.tblout').name
        # Initialize domtblout temporary file path
        domtblout_path = NamedTemporaryFile(delete=False, dir=os.path.dirname(hmm_path), suffix='.domtblout').name

        # Safely run hmmsearch
        try:
            # Run hmmsearch against given chunk
            hmm_search.run(
                # Path to searched HMM model/library file
                hmm_path=hmm_path,
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

        # Catch exception
        except CalledProcessError:
            # Remove unzipped temporary file
            if chunk_gzip:
                os.remove(chunk_path)
            # Raise the error interrupting pipeline
            raise

        # Remove unzipped temporary file
        if chunk_gzip:
            os.remove(chunk_path)

        # Return domtblout path
        return tblout_path, domtblout_path

    # Search HMM against a target dataset
    def search_hmm_dataset(self, hmm_path, target_dataset, e_value, z_score, n_cores, client, verbose=False):
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
        target_dataset (list)   List of target dataset chunks
        e_value (float)         E-value to define a match significant
        z_score (int)           Z-score defining whole dataset length
        n_cores (int)           Number of cores available per job
        client (Client)         Client used to spawn jobs on cluster
        verbose (bool)          Whether to print verbose log or not

        Return
        (list)                  List of TBLOUT tabular output files
        (list)                  List of DOMTBLOUT tabular output files

        Raise
        (CalledProcessError)    In case hmmsearch could not run
        """
        # Run futures safely
        try:

            # Initialize futures
            futures = list()
            # Loop through each target dataset
            for chunk in target_dataset:
                # Get chunk path
                chunk_path = chunk.path
                # Submit search function
                futures.append(client.submit(
                    func=self.search_hmm_chunk,
                    hmm_path=hmm_path,
                    chunk_path=chunk_path,
                    e_value=e_value,
                    z_score=z_score,
                    n_cores=n_cores,
                    hmm_search=self.hmm_search
                ))

            # Retrieve results (tuples containing paths to tabluar output files)
            results = client.gather(futures)
            # Define output lists
            tblout_list = [results[i][0] for i in range(len(results))]
            domtblout_list = [results[i][1] for i in range(len(results))]
            # Return retrieved TBLOUT and DOMTBLOUT lists
            return tblout_list, domtblout_list

        except CalledProcessError as err:
            # Verbose log
            if verbose:
                # Print error STDERR
                print(' '.join([
                    'Could not run distributed hmmsearch:',
                    'returned code {}:\n'.format(err.returncode),
                    err.stdout.trim()
                ]))
            # Raise error interrupting pipeline
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
        client(Client)              Dask Client used to submit distributed jobs

        Return
        (dict)      Dictionary mapping cluster name to sequences bit-score thresholds
        (dict)      Dictionary mapping cluster name to domains bit-score thresholds
        (dict)      Dictionary mapping cluster name to sequences hits
        (dict)      Dictionary mapping cluster name to domains hits
        """
        # # Initialize results list
        # results = list()
        # # Initialize futures
        # futures = list()
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
            # # Distribute scores retrieval function
            # futures.append(client.submit(
            #     func=self.parse_search_scores,
            #     tblout_type='sequences',
            #     tblout_path=tblout_path,
            #     e_value=e_value,
            #     min_bits=25.0,
            #     max_bits=999999.99
            # ))
        # # Gather results
        # results = client.gather(futures)
        # # Merge all retrieved scores
        # sequence_scores = merge_scores(*results)

        # # Reinitialize futures
        # futures = list()
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
            # # Distribute hits retrieval function
            # futures.append(client.submit(
            #     func=self.parse_search_hits,
            #     tblout_type='sequences',
            #     tblout_path=tblout_path,
            #     scores=sequence_scores
            # ))
        # # Gather results (list of dictionaries mapping clusters to sequences)
        # results = client.gather(futures)
        # # Merge all retrieved dictionaries
        # sequence_hits = merge_hits(*results)

        # # Reinitialize futures
        # futures = list()
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
            # # Distribute scores retrieval function
            # futures.append(client.submit(
            #     func=self.parse_search_scores,
            #     tblout_type='domains',
            #     tblout_path=domtblout_path,
            #     e_value=e_value,
            #     min_bits=25.0,
            #     max_bits=999999.99
            # ))
        # # Gather results
        # results = client.gather(futures)
        # # Merge all retrieved scores
        # domain_scores = merge_scores(*results)

        # # Reinitialize futures
        # futures = list()
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
            # Keep only domains which have also a hit in sequences
            hits = {
                qname: hits[qname] for qname
                in (set(hits.keys()) & set(sequence_hits.keys()))
            }
            # Update domain hits
            domain_hits = merge_hits(domain_hits, hits)
            # # Distribute hits retrieval function
            # futures.append(client.submit(
            #     func=self.parse_search_hits,
            #     tblout_type='domains',
            #     tblout_path=domtblout_path,
            #     scores=domain_scores
            # ))
        # # Gather results (list of dictionaries mapping clusters to sequences)
        # results = client.gather(futures)
        # # Merge all retrieved dictionaries
        # domain_hits = merge_hits(*results)

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
        # Define cluster iterator
        clusters_iter = glob(os.path.join(clusters_path, 'MGYP*'))
        # Loop through each cluster in clusters path
        cluster_names = {
            os.path.basename(cluster_path)
            for cluster_path
            in clusters_iter
        }

        # # Initialize list of available HMM models
        # hmm_models = list()
        # # Initialize length of longest model
        # hmm_length = 0
        # # Retrieve available HMM models
        # for cluster_path in clusters_iter:
        #     # Get cluster name
        #     cluster_name = os.path.basename(cluster_path)
        #     # Retrieve HMM path
        #     hmm_path = os.path.join(cluster_path, 'HMM')
        #     # Ensure file exists
        #     if not os.path.isfile(hmm_path):
        #         # Eventually raise new exceptio
        #         raise FileNotFoundError(' '.join([
        #             'Error: could not find HMM model',
        #             'for cluster {:s}'.format(cluster_name),
        #             'at {:s}'.format(hmm_path)
        #         ]))
        #     # Load HMM model
        #     hmm_model = HMM.from_file(hmm_path)
        #     # Update longest
        #     hmm_length = max(hmm_length, hmm_model.length)
        #     # Store HMM model
        #     hmm_models.append(hmm_model)

        # Make HMM library
        hmmlib_path, hmmlib_len, cluster_names = self.make_hmm_library(
            cluster_names=cluster_names,
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Get maximum number of cores
        max_cores = int(self.cluster_kwargs.get('cores', 1))
        # Get maximum available memory per job
        max_memory = as_bytes(self.cluster_kwargs.get('memory', '0 Gb'))

        # Verbose log
        if verbose:
            print('Maximum resources allowed per job are', end=' ')
            print('{:d} cores and'.format(max_cores), end=' ')
            print('{:d} Gb of memory'.format(math.floor(max_memory / 1e09)))

        # Get minimum required memory per CPU
        req_memory = HMMSearch.get_memory(
            model_len=hmmlib_len,
            longest_seq=uniprot_shape[1]
        )
        # Define maximum number of CPUs (cores) allocable to run hmmsearch
        num_cores = HMMSearch.get_cpus(
            model_len=hmmlib_len,
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

        # Verbose log
        if verbose:
            print('Making HMM search against UniProt dataset', end=' ')
            print('on {:d} cores with'.format(num_cores), end=' ')
            print('{:d} Gb of memory each'.format(math.floor(max_memory / num_cores / 1e09)))

        # Instantiate new cluster with given parameters
        cluster, client = self.get_client(dict(cores=num_cores))
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Make search against UniProt (retrieve tabular output files)
        time_run, (tblout_paths, domtblout_paths) = benchmark(
            fn=self.search_hmm_dataset,
            hmm_path=hmmlib_path,
            target_dataset=self.uniprot,
            e_value=1000,
            z_score=uniprot_shape[0],
            n_cores=num_cores,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print(' '.join([
                'HMM models search against UniProt done',
                'in {:.0f} seconds'.format(time_run)
            ]))

        # # DEBUG
        # print('DEBUG')
        # print(tblout_paths[:10])
        # print(domtblout_paths[:10])

        # Close client
        self.close_client(cluster, client)

        # Open new client with less resources
        cluster, client = self.get_client(dict(cores=1))
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Parse tabular output files: retrieve significant scores and hits
        time_run, (sequence_scores, sequence_hits, domain_scores, domain_hits) = benchmark(
            fn=self.parse_hmm_hits,
            tblout_paths=tblout_paths,
            domtblout_paths=domtblout_paths,
            e_value=self.uniprot_e_value,
            client=client
        )

        # # DEBUG
        # print('Sequence scores:', sequence_scores)
        # print('Domain scores:', domain_scores)
        # print('Sequence_hits:', sequence_hits)
        # print('Domain hits:', domain_hits)

        # Go through every temporary output files
        for temp_path in [*tblout_paths, *domtblout_paths]:
            # Remove temporary file
            os.remove(temp_path)

        # Define cluster names in sequences and domain hits
        intersect_uniprot = set(sequence_hits.keys()) | set(domain_hits.keys())
        # Udpate cluster names set
        cluster_names = cluster_names - intersect_uniprot
        # Loop through each cluster intersecting uniprot
        for cluster_name in intersect_uniprot:
            # Define current path
            curr_path = os.path.join(clusters_path, cluster_name)
            # Define target path
            move_path = os.path.join(uniprot_path, cluster_name)
            # Move cluster to uniprot directory
            shutil.move(curr_path, move_path)

        # Verbose log
        if verbose:
            print(' '.join([
                'HMM models hits against UniProt retrieved',
                'in {:.0f} seconds:'.format(time_run),
                '{:d} clusters have been kept and'.format(len(cluster_names)),
                '{:d} removed'.format(len(intersect_uniprot)),
                'due to UniProt intersections'
            ]))

        # Close client
        self.close_client(cluster, client)

    # Search multiple HMM against MGnifam
    def search_hmm_mgnifam(self, clusters_path, mgnifam_shape, min_jobs=1, max_jobs=100, verbose=False):
        """ Serach multiple HMM models against MGnifam

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

        # Define cluster iterator
        clusters_iter = glob(os.path.join(clusters_path, 'MGYP*'))
        # Loop through each cluster in clusters path
        cluster_names = {
            os.path.basename(cluster_path)
            for cluster_path
            in clusters_iter
        }

        # Make HMM library
        hmmlib_path, hmmlib_len, cluster_names = self.make_hmm_library(
            cluster_names=cluster_names,
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Get maximum number of cores
        max_cores = int(self.cluster_kwargs.get('cores', 1))
        # Get maximum available memory per job
        max_memory = as_bytes(self.cluster_kwargs.get('memory', '0 Gb'))
        # Verbose log
        if verbose:
            print('Maximum resources allowed per job are', end=' ')
            print('{:d} cores and'.format(max_cores), end=' ')
            print('{:d} Gb of memory'.format(math.floor(max_memory / 1e09)))

        # Get minimum required memory per CPU
        req_memory = HMMSearch.get_memory(
            model_len=hmmlib_len,
            longest_seq=mgnifam_shape[1]
        )
        # Define maximum number of CPUs (cores) allocable to run hmmsearch
        num_cores = HMMSearch.get_cpus(
            model_len=hmmlib_len,
            max_memory=max_memory,
            max_cpus=max_cores,
            longest_seq=mgnifam_shape[1]
        )
        # Check if memory is sufficient
        if not num_cores:
            # Raise new memory exception
            raise MemoryError(' '.join([
                'Unable to start HMM search against MGnifam dataset:',
                'minimum required memory is {:d} GB'.format(math.ceil(req_memory / 1e09)),
                'while maximum allocable is {:d} GB'.format(math.floor(max_memory / 1e09))
            ]))
        # Verbose log
        if verbose:
            print('Making HMM search against Mgnifam dataset:', end=' ')
            print('on {:d} cores with'.format(num_cores), end=' ')
            print('{:d} Gb of memory each'.format(math.floor(max_memory / 1e09)))

        # Instantiate new cluster with given parameters
        cluster, client = self.get_client(dict(cores=num_cores))
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Make search against UniProt (retrieve tabular output files)
        time_run, (tblout_paths, domtblout_paths) = benchmark(
            fn=self.search_hmm_dataset,
            hmm_path=hmmlib_path,
            target_dataset=self.mgnifam,
            e_value=1000,
            z_score=mgnifam_shape[0],
            n_cores=num_cores,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print(' '.join([
                'Search against MGnifam done',
                'in {:.0f} seconds'.format(time_run)
            ]))

        # Close client
        self.close_client(cluster, client)

        # Open new client with less resources
        cluster, client = self.get_client(dict(cores=1))
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Parse tabular output files: retrieve significant scores and hits
        time_run, (sequence_scores, sequence_hits, domain_scores, domain_hits) = benchmark(
            fn=self.parse_hmm_hits,
            tblout_paths=tblout_paths,
            domtblout_paths=domtblout_paths,
            e_value=self.mgnifam_e_value,
            client=client
        )

        # Remove temporary output files
        for temp_path in [*tblout_paths, *domtblout_paths]:
            # Remove tabular output file
            os.remove(temp_path)

        # # DEBUG Show domain hits
        # for i, hit in enumerate(domain_hits):
        #     # Print current hit
        #     print('{:d}-th domain hit:'.format(i+1), hit)
        #     # Early stopping
        #     if i >= 10:
        #         break


        # Initialize sequences to retrieve
        sequences_acc = set()
        # Loop through each sequence accession in domain hits
        for cluster_name in domain_hits.keys():
            # Loop through each sequence accession
            for sequence_acc in domain_hits[cluster_name]:
                # Store sequence accession
                sequences_acc.add(sequence_acc)

        # # DEBUG Show sequences accession
        # for i, sequence_acc in enumerate(sequences_acc):
        #     # Print current hit
        #     print('{:d}-th sequence accession:'.format(i+1), sequence_acc)
        #     # Early stopping
        #     if i >= 10:
        #         break

        # Fill mgnifam dictionary
        time_run, (fasta_sequences, num_sequences) = benchmark(
            fn=self.search_fasta_sequences,
            sequences_acc=sequences_acc,
            fasta_dataset=self.mgnifam,
            client=client,
            verbose=verbose
        )

        # # DEBUG
        # # Show number of sequences checked
        # print('Number of retrieved sequences:', num_sequences)
        # # Store fasta sequences file
        # with open(os.path.join(clusters_path, 'fasta_sequences.json'), 'w') as f:
        #     # Write to file
        #     json.dump(fasta_sequences, f, indent=2)

        # Align cluster sequences to HMM using hmmalign
        self.make_hmm_alignments(
            clusters_path=clusters_path,
            fasta_sequences=fasta_sequences,
            cluster_members=domain_hits,
            client=client
        )

        # Close client
        self.close_client(cluster, client)

    # Make HMM alignment
    def make_hmm_alignments(self, clusters_path, fasta_sequences, cluster_members, client, verbose=False):
        """ Align all HMM models to sequneces

        Args
        Return
        Raise
        """
        # Initialize futures
        futures = list()
        # Loop through each cluster path
        for cluster_path in glob(os.path.join(clusters_path, 'MGYP*')):
            # Define cluster name
            cluster_name = os.path.basename(cluster_path)
            # Define fet of sequence accession from cluster members
            cluster_acc = cluster_members.get(cluster_name, set())

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
            futures.append(client.submit(
                func=self.align_hmm_cluster,
                cluster_path=cluster_path,
                fasta_path=fasta_path,
                # fasta_sequences=fasta_sequences,
                # sequences_acc=cluster_acc,
                hmm_align=self.hmm_align
            ))
        # Gather all futures
        client.gather(futures)

    @staticmethod
    def align_hmm_cluster(cluster_path, fasta_path, hmm_align):
        """ Align HMM to cluster sequneces

        Raise
        (KeyError)          In case fasta entry can not be retrieved
        """

        # Define cluster name
        cluster_name = os.path.basename(cluster_path)

        # Define input fasta file path
        fasta_path = os.path.join('ALIGN.fasta')

        # Define output alignment path (.sto stockholm)
        sto_path = os.path.join(cluster_path, 'ALIGN.sto')
        # Define output  alignment path (.aln)
        aln_path = os.path.join(cluster_path, 'ALIGN')

        # Define input HMM model path
        hmm_path = os.path.join(cluster_path, 'HMM')
        # Check if HMM model path is valid
        if not os.path.isfile(hmm_path):
            # Raise new exception
            raise FileNotFoundError(' '.join([
                'Unable to run hmmalign script:',
                'HMM model for cluster {:s}'.format(cluster_name),
                'not found at {:s}'.format(hmm_path)
            ]))

        # # Define fasta path
        # fasta_path = NamedTemporaryFile(delete=False, suffix='.fa').name
        # # Fill fasta file with sequences
        # with open(fasta_path, 'w') as fasta_file:
        #     # Loop through each cluster members
        #     for sequence_acc in sequences_acc:
        #         # Retrieve fasta entry
        #         fasta_entry = fasta_sequences.get(sequence_acc)
        #         # Case entry is empty
        #         if not fasta_entry:
        #             # Raise exception
        #             raise KeyError(' '.join([
        #                 'Error: could not retrieve fasta sequence',
        #                 '{:s}'.format(sequence_acc),
        #                 'for cluster {:s}'.format(cluster_name)
        #             ]))
        #         # Add entry to file
        #         fasta_file.write(fasta_entry + '\n')

        # Run hmmalign script
        hmm_align.run(
            # Set path to input HMM model
            hmm_path=hmm_path,
            # Set path to input FASTA file
            fasta_path=fasta_path,
            # Set path to output alignment file
            out_path=sto_path
        )
        # # Remove temporary fasta file
        # os.remove(fasta_path)

        # # TODO Load retrieved alignment file
        # align = MSA.from_sto(sto_path)
        # # TODO Parse loaded alignment file
        # align.to_aln(aln_path)


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


# Test
if __name__ == '__main__':

    # Define root path
    ROOT_PATH = os.path.dirname(__file__) + '/../..'

    # Initialize scores
    scores = list()
    # Define TBLOUT iterator
    tblout_iter = glob(ROOT_PATH + '/tmp/release2/batches/0/MGYP*/*.domtblout')
    # # DEBUG
    # print(tblout_iter)
    # Go through all sample tblout paths
    for i, tblout_path in enumerate(tblout_iter):
        # Retrieve scores
        scores.append(Build.parse_search_scores(
            tblout_path=tblout_path,
            tblout_type='domains',
            e_value=0.001,
            min_bits=25.0,
            max_bits=999999.99
        ))
        # # Print scores
        # print('{:d}-th scores:\n'.format(i+1), scores[-1])
        # print()
    # Merge scores together
    scores = merge_scores(*scores)

    # Show results
    print('Scores:')
    for cluster_name, cluster_score in scores.items():
        print(cluster_name, cluster_score)
    print()

    # Initialize hits
    hits = list()
    # Go through all sample tblout paths
    for tblout_path in tblout_iter:
        # Retrieve hits
        hits.append(Build.parse_search_hits(
            tblout_path=tblout_path,
            tblout_type='domains',
            scores=scores
        ))
        # # Print latest retrieved hits
        # print('{:d}-th hits\n'.format(), hits[-1])
        # print()
    # Merge hits together
    hits = merge_hits(*hits)

    # Show results
    print('Hits:')
    for cluster_name, cluster_members in hits.items():
        print(cluster_name, cluster_members)
    print()
