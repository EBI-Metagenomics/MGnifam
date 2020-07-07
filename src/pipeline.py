# Dependencies
from dask.distributed import Client
from tqdm import tqdm
from glob import glob
import numpy as np
import time
import json
import sys
import os
import re

# Plotting libraries
import matplotlib
import matplotlib.pyplot as plt
# Set plotting non-interactive
matplotlib.use('Agg')

# Custom dependencies
import src.dataset as ds
import src.disorder as disorder
import src.msa as msa
import src.hmm as hmm


class Pipeline(object):

    # Constructor
    def __init__(
        self, in_path, out_path, dask_client,
        batch_size=100, num_clusters=1000,
        clusters_path='data/clusters/chunk*.tsv.gz',
        mgnify_path='data/mgnify/chunk*.fa.gz'
    ):

        # Check input path type
        self.in_path = self.check_path(in_path)

        # Case output directory does not exist
        if not os.path.isdir(out_path):
            # Make output directory
            os.mkdir(out_path)
        # Set output directory and ensure format
        self.out_path = os.path.dirname(out_path)

        # Save batch parameters
        self.batch_size = batch_size
        self.num_clusters = num_clusters

        # Check path to clusters dataset
        self.clusters_path = self.check_path(clusters_path)
        # Check path to mgnify dataset
        self.mgnify_path = self.check_path(mgnify_path)

        # Set client (Client(LocalCluster) by default)
        self.dask_client = dask_client

    # Wrapper for run method
    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    # Run abstract method
    def run(self, *args, **kwargs):
        raise NotImplementedError

    # Check list of paths / single path
    @staticmethod
    def check_path(path):
        # Check input path type: string
        if isinstance(path, str):
            # Use glob to find files referenced by unix string
            return glob(path)
        # Case input path is already a list
        elif isinstance(in_path, list):
            # Return given list
            return path
        # Otherwise, raise an error
        else:
            raise ValueError('Given path is not valid')

    # Update a log dictionary, eventually writing it out to a given path
    @staticmethod
    def update_log(log_old, log_new, log_path=None):
        # Update log dictionary
        for k, v in log_new.items():
            log_old[k] = v
        # Case log path is set
        if isinstance(log_path, str) and log_path:
            # Open file buffer
            with open(log_path, 'w') as lf:
                # Dump dictionary
                json.dump(log_old, lf, indent=4)
        # return new dictionary
        return log_old


class Release(Pipeline):

    def __init__(self, *args, **kwargs):
        raise NotImplementedError


class Seed(Pipeline):

    # Constructor
    def __init__(
        self, clusters_path, mgnify_path, dask_client,
        mobidb_env=os.environ.copy(), mobidb_cmd=['mobidb_lite.py'],
        muscle_env=os.environ.copy(), muscle_cmd=['muscle'],
        hmm_build_env=os.environ.copy(), hmm_build_cmd=['hmmbuild'],
        hmm_search_env=os.environ.copy(), hmm_search_cmd=['hmmsearch']
    ):
        # Retrieve clusters dataset
        clusters_path = self.check_path(clusters_path)
        self.ds_clusters = ds.Cluster.from_list(clusters_path)
        # Retrieve MGnify dataset
        mgnify_path = self.check_path(mgnify_path)
        self.ds_mgnify = ds.MGnify.from_list(mgnify_path)
        # Initialize a new MobiDB Lite instance
        self.mobidb = disorder.MobiDBLite(cmd=mobidb_cmd, env=mobidb_env)
        # Initialize a new Muscle instance
        self.muscle = msa.Muscle(cmd=muscle_cmd, env=muscle_env)
        # Initialize a new hmmbuild instance
        self.hmm_build = hmm.HMMBuild()
        # Initialize a new hmmsearch instance
        self.hmm_search = hmm.HMMSearch()
        # Set client (Client(LocalCluster) by default)
        self.dask_client = dask_client

    # Run Seed alignment (for a single cluster)
    def run(self, cluster_names, cluster_dir, verbose=False, log=False):

        # Initialize log
        log = {
            # Path to output clusters directory
            'cluster_path': '',
            # Number of input cluster names
            'num_clusters': 0,
            # Number of sequences
            'num_sequences': 0,
            # Number of sequences per cluster
            'cluster_sizes': 0,
            # Compositional bias dict(cluster name: compositional bias)
            'comp_bias': dict(),
            # Define total running time (seconds)
            'time_took': 0.0
        }
        # Initialize timers
        time_beg, time_end = time.time(), None
        # Define log path (only if required)
        log_path = os.path.join(cluster_dir, 'log.json') if log else ''

        # Case output directory does not exits
        if not os.path.exists(cluster_dir):
            # Make directory
            os.mkdir(cluster_dir)
        # Update log
        self.update_log(log_path=log_path, log_old=log, log_new={
            'cluster_path': cluster_dir
        })

        # Retrieve cluster members from clusters .tsv dataset
        cluster_members = self.get_cluster_members(
            cluster_names=cluster_names,
            verbose=verbose
        )
        # Update log
        self.update_log(log_path=log_path, log_old=log, log_new={
            'num_clusters': len(cluster_members),
            'cluster_sizes': np.mean([
                len(cluster)  # Take size of the current cluster
                for cluster  # Get each cluster in cluster members
                in cluster_members.values()  # Get cluster sequences
            ])
        })

        # Check cluster members raises error
        self.check_cluster_members(cluster_members)

        # Define sequences accession numbers list
        sequences_acc = set([
            # Get i-th sequence accession in cluster named n
            cluster_members[n][i]
            # For each cluster in cluster members
            for n in cluster_members
            # For each component of that cluster
            for i in range(len(cluster_members[n]))
        ])
        # Update log
        self.update_log(log_path=log_path, log_old=log, log_new={
            'num_sequences': len(sequences_acc)
        })

        # Retrieve fasta sequences from MGnify .fa dataset
        fasta_sequences = self.get_fasta_sequences(
            sequences_acc=sequences_acc,
            verbose=verbose
        )

        # Check for mismatching clusters wrt retrieved sequences
        self.check_cluster_sequences(cluster_members, fasta_sequences)

        # Compute compositional bias
        comp_bias = self.compute_comp_bias(
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            verbose=verbose
        )
        # Update log
        self.update_log(log_path=log_path, log_old=log, log_new={
            'comp_bias': np.mean([*comp_bias.values()])
        })

        # Plot compositional bias
        self.plot_comp_bias(out_path=os.path.join(cluster_dir, 'comp_bias.png'))

        # TODO Remove clusters with compositional bias too high

        # Loop through each cluster name to make clusters directories
        for cluster_name in cluster_members:
            # Define cluster subdir
            cluster_subdir = os.path.join(cluster_dir, cluster_name)
            # Check if a subdirectory for current cluster exists
            if not os.path.exists(cluster_subdir):
                # Make sub directory
                os.mkdir(cluster_subdir)

        # Make seed alignments
        self.make_sequences_aln(
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            cluster_dir=cluster_dir,
            verbose=verbose
        )

        # Make HMM
        self.build_hmms(
            cluster_members=cluster_members,
            cluster_dir=cluster_dir,
            verbose=verbose
        )

        # TODO Search HMM against MGnify

        # Update timers
        time_end = time.time()
        # Save run duration (seconds)
        self.update_log(log_path=log_path, log_old=log, log_new={
            'time_took': float(time_end - time_beg)
        })

    # Retrieve cluster members for a batch of clusters
    def get_cluster_members(self, cluster_names, verbose=False):
        # Initialize output dict(cluster name: sequence acc)
        cluster_members = dict()
        # Intialize timers
        time_beg, time_end = time.time(), None

        # Verbose log
        if verbose:
            # Get a summary of cluster names
            print('Retrieving {:d} clusters:\n{}'.format(
                len(cluster_names),  # Number of clusters
                ', '.join(cluster_names)  # Verbose list of cluster names
            ))

        # Submit jobs for searching cluster members
        futures = self.dask_client.map(
            # Search for cluster members in current chunk
            lambda ds: ds.search(cluster_names),
            # List of cluster datasets
            self.ds_clusters
        )
        # Retrieve found clusters (list of lists)
        results = self.dask_client.gather(futures)
        # Loop through resulting clusters dictionaries
        for i in range(len(results)):
            # Loop through every cluster name in result dictionary
            for cluster_name in results[i].keys():
                # Initialize entry for current cluster
                cluster_members.setdefault(cluster_name, [])
                # Update cluster members list
                cluster_members[cluster_name] += results[i][cluster_name]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Get execution time summary
            print('Took {:.0f} seconds to search members of {:d} clusters'.format(
                time_end - time_beg,  # Time required to distributed search
                len(cluster_members.keys())  # Number of clusters
            ))

        # Return cluster members dict
        return cluster_members

    # Check that all clusters have at least one sequence accession
    def check_cluster_members(self, cluster_members):
        # Check for clusters with no members (eventually raise error)
        empty_clusters = [n for n, m in cluster_members.items() if not len(m)]
        # Check if at least one cluster has no member
        if empty_clusters:
            # Raise new error
            raise KeyError('No sequences found for {:d} clusters: {:s}'.format(
                len(empty_clusters),  # Number of clusters with no members
                ', '.join(empty_clusters)  # Verbose string of empty clusters
            ))

    # Retrieve cluster sequences for given sequences accession numbers
    def get_fasta_sequences(self, sequences_acc, verbose=False):
        # Initialize output dict(sequence acc: fasta entry)
        fasta_sequences = dict()
        # Initialize timers
        time_beg, time_end = time.time(), None

        # Submit jobs for searching cluster members sequences
        futures = self.dask_client.map(
            # Search for cluster members in current chunk
            lambda ds: ds.search(sequences_acc),
            # List of cluster datasets
            self.ds_mgnify
        )
        # Retrieve found sequences (list of dict)
        results = self.dask_client.gather(futures)
        # Loop through resulting sequences dictionaries
        for i in range(len(results)):
            # Loop through retrieved sequence accession numbers
            for acc in results[i]:
                # Update cluster sequences dictionary
                fasta_sequences[acc] = results[i][acc]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Log time taken for fasta sequences distributed search
            print('Took {:.0f} seconds to search for {:d} sequences'.format(
                time_end - time_beg,  # Time required to distributed search
                len(fasta_sequences)  # Number of sequences
            ))

        # Return fasta sequences
        return fasta_sequences

    # Check that all the sequences have been retrieved
    def check_cluster_sequences(self, cluster_members, fasta_sequences):
        # Initialize list of clusters with mismatching sequences
        mismatch_clusters = list()
        # Get sequence accession numbers in retrieved fasta sequences
        sequences_acc = set(fasta_sequences.keys())
        # Fill mismatching clusters
        for cluster_name in cluster_members.keys():
            # Define current cluster's seuqence accessionnumbers
            cluster_acc = set(cluster_members[cluster_name])
            # Case some sequence has not been retrieved
            if cluster_acc - sequences_acc:
                # Add cluster name to list of mismatching ones
                mismatch_clusters.append(cluster_name)
        # Case there is at least one mismatching cluster
        if mismatch_clusters:
            # Raise new error
            raise KeyError(
                'Could not retrieve all sequences for {:d} clusters: {}'.format(
                    len(mismatch_clusters),  # Number of mismatching clusters
                    ', '.join(mismatch_clusters)  # Mismatching clusters
                )
            )

    @staticmethod
    def compute_comp_bias_(sequences, mobidb, threshold=1, inclusive=True):
        """Compute compositional bias of a single cluster
        First, calls MobiDB Lite to predict disordered regions
        Then, applies a threshold over minimum number of predictions required
        to define a reidue actually disordered.
        Lastly, computes compositional bias over the thesholded predictions.

        Args
        sequences (dict(str: str))      Dictionary linking each cluster member
                                        sequence accession number to the
                                        associated fasta sequence
        mobidb (disorder.MobiDBLite)    Instance of MobiDB Lite
        threshold (int/str)             Threshold on sum of disorder predictions
                                        in order to define a single residue
                                        actually disordered
        inclusive (bool)                Wether the threshold must be
                                        inclusive (i.e. #pred >= threshold) or
                                        exclusive (i.e. #pred > threshold)

        Return
        (float)                         Compositional bias
        """
        # Make disordered regions predictions
        pred = mobidb.run(sequences=sequences)
        # Apply threshold
        pred = disorder.compute_threshold(
            sequences=pred,
            threshold=threshold,
            inclusive=inclusive
        )
        # Compute compositional bias
        return disorder.compute_comp_bias(pred)

    # Compute compositional bias for every sequence
    def compute_comp_bias(self, cluster_members, fasta_sequences, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None

        # Initialize output dict(cluster name: compositional bias)
        comp_bias = dict()

        # Get cluster names
        cluster_names = [*cluster_members.keys()]

        # Initialize futures list
        futures = list()
        # Loop through every cluster name
        for cluster_name in cluster_names:
            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = dict()
            # Loop through each sequence accession in cluster
            for acc in cluster_members[cluster_name]:
                # Save fasta sequence
                cluster_sequences[acc] = fasta_sequences[acc]
            # Run distributed compositional bias
            futures.append(self.dask_client.submit(
                self.compute_comp_bias_,
                sequences=cluster_sequences,
                mobidb=self.mobidb,
                threshold=1,
                inclusive=True
            ))
        # Get results
        results = self.dask_client.gather(futures)
        # Associate cluster name to each compositional bias
        for i in range(len(cluster_names)):
            # Save compositional bias for current cluster
            comp_bias[cluster_names[i]] = results[i]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Get execution time summary
            print('Took {:.0f} seconds to compute disordered regions for {:d} fasta sequences'.format(
                time_end - time_beg,  # Time required to distributed search
                len(comp_bias)  # Number of fasta sequences
            ))

        # Return dict(cluster name: compositional bias)
        return comp_bias

    # Make a plot of compositional bias
    def plot_comp_bias(self, comp_bias, out_path):
        # Plot compositional bias distribution
        fig, ax = plt.subplots(figsize=(10, 5))
        # Make plot
        ax.set_title('Compositional Bias distribution')
        ax.set_xlabel('Compositional Bias')
        ax.set_ylabel('Number of clusters')
        ax.hist(comp_bias.values(), density=False, bins=100)
        ax.set_xlim(left=0.0, right=1.0)
        # Save plot
        plt.savefig(out_path)
        # Close plot
        plt.close()

    @staticmethod
    def make_sequences_aln_(sequences, muscle, aln_path):
        """Aligns a list of fasta sequences
        Takes as input a list of sequences and aligns them, then store the
        alignment as a .aln file.

        Args
        sequences (dict(str: str))
        muscle (msa.Muscle)         Multiple Sequence Alignment (MSA) algorithm
                                    used for generating output .aln file
        path (str)                  Path to output multiple sequence alignment
                                    file
        """
        # Make multiple sequence alignment
        aln = muscle.run(sequences=sequences, acc_regex='>(\S+)', verbose=False)
        # Store MSA to disk
        aln.to_aln(out_path=aln_path)
        # Free memory
        del aln

    # Make multiple sequence alignment
    def make_sequences_aln(self, cluster_members, fasta_sequences, cluster_dir, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None
        # Verbose log
        if verbose:
            print('Making multiple sequence alignments for {:d} clusters...'.format(
                len(cluster_members)
            ))

        # Get cluster names
        cluster_names = [*cluster_members.keys()]

        # Initialize futures list
        futures = list()
        # Loop through every cluster name
        for cluster_name in cluster_names:
            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = list()
            # Loop through each sequence accession in cluster
            for acc in cluster_members[cluster_name]:
                # Save fasta sequence
                cluster_sequences.append(fasta_sequences[acc])
            # Run distributed compositional bias
            futures.append(self.dask_client.submit(
                self.make_sequences_aln_,
                sequences=cluster_sequences,
                muscle=self.muscle,
                aln_path=os.path.join(cluster_dir, cluster_name, 'SEED')
            ))
        # Get results
        self.dask_client.gather(futures)

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            print('Took {:.0f} seconds to make {:d} multiple sequence alignments'.format(
                time_end - time_beg,
                len(cluster_members)
            ))

    # Create Hidden Markov Models (HMMs)
    def build_hmms(self, cluster_members, cluster_dir, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None
        # Verbose log
        if verbose:
            print('Building Hidden Markov Models (HMMs) for {:d} clusters...'.format(
                len(cluster_members)
            ))

        # Intialize futures container
        futures = list()
        # Loop through each cluster name
        for cluster_name in cluster_members:
            # Define path to current cluster
            cluster_path = os.path.join(cluster_dir, cluster_name)
            # Define path to current cluster seed alignment
            seed_path = os.path.join(cluster_path, 'SEED')
            # Define path to output Hidden Markov Model (HMM)
            hmm_path = os.path.join(cluster_path, 'HMM')
            # Run HMM build (distributed)
            futures.append(self.dask_client.submit(
                # Function to be submitted
                self.hmm_build.run,
                # Parameters to be feed to submitted function
                msa_path=seed_path,  # Set path to input seed alignment
                out_path=hmm_path,  # Set path to output file
                name=cluster_name  # Set HMM name as current cluster name
            ))

        # Retrieve results
        self.dask_client.gather(futures)

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            print('Took {:.0f} seconds to build {:d} Hidden Markov Models (HMMs)'.format(
                time_end - time_beg,
                len(cluster_members)
            ))

    @staticmethod
    def search_hmms_(hmm_search, hmm_path, target_ds):

        # Initialize a temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Search given hmm aganist given dataset using hmmsearch
        hmm_search.run(
            # Path to searched HMM file
            hmm_path=hmm_path,
            # Path to sequences target dataset
            ds_path=target_ds.path,
            # Path to per-domain output table
            domtblout_path=out_file.name
        )

        # Initialize domtblout
        domtblout = list()
        # Open output temporary file
        with open(out_file.name) as of:
            # Save domtblout as list of rows
            for row in hmm_search.iter_domtblout(of):
                # Store current row
                domtblout.append(row)

        # Remove temporary output file
        os.remove(out_file.name)

        # Return per-domain table
        return domtblout

    # Search Hidden Markov Models (HMMs) against a dataset
    def search_hmms(self, cluster_members, cluster_dir, target_ds, ds_size, verbose=False):
        """Search multiple HMMs against given dataset

        Args
        cluster_members (iterable)  Contains names of clusters whose HMM must be
                                    searched against given dataset
        cluster_dir (str)           Root directory of clusters sub-directories
        target_ds (list)            List of dataset chunks against which the
                                    clusters Hidden Markov Models must be run
                                    in a distributed way
        ds_size (int)               Number of entries in input dataset, in total
        verbose (bool)              Wether to return or not verbose output

        Return
        (list)                      Concatenated domtblout table, containing
                                    all matches between HMMs and sequences
        """

        # Initialize a single HMM temporary file containing all built HMMs
        hmm_temp = tempfile.NamedTemporaryFile(delete=False)
        # Open write buffer to full HMM file (append)
        hmm_file = open(hmm_temp.name, 'a')

        # Loop through each cluster names
        for cluster_name in cluster_members:

            # Define full path to current HMM file
            curr_path = os.path.join(custer_dir, cluster_name, 'HMM')
            # Retrieve current hmm file
            with open(curr_path, r) as curr_file:
                # Loop thriugh every line of current HMM file
                for line in curr_file:
                    # Append line to full hmm file
                    hmm_file.write(line)
        # Close write buffer to full HMM file
        hmm_file.close()

        # Initialize futures container
        futures = list()
        # Loop through each dataset in <against_ds> list
        for i in range(len(target_ds)):
            # Distribute HMM search against target dataset
            self.dask_client.submit(
                # Function to submit
                self.search_hmms_,
                # Function parameters
                hmm_search=self.hmm_search,
                hmm_path=hmm_file.name,
                target_ds=target_ds[i]
            )
        # Retrieve results
        results = self.dask_client.gather(futures)
        # Concatenate all results in a single table
        domtblout = [
            # Get j-th domtblout row of i-th dataset chunk
            results[i][j]
            # Loop through each dataset chunk result
            for i in range(len(results))
            # Loop through each row in i-th dataset chunk result
            for j in range(len(results[i]))
        ]
        # Delete results (free space)
        del results

        # Delete temporary file
        os.remove(hmm_temp.name)

    # Wrapper for searching given HMMs against Mgnify dataset
    def search_hmms_mgnify(self, *args, **kwargs):
        # Call generic function
        return self.search_hmms(*args, **kwargs, target_ds=self.ds_mgnify)

    # # TODO Wrapper for searching given HMMs against UniProt dataset
    # def search_hmms_mgnify(self, *args, **kwargs):
    #     # Call generic function
    #     return self.search_hmms(*args, **kwargs, target_ds=self.ds_uniprot)
