# Dependencies
from dask.distributed import Client
from tqdm import tqdm
from glob import glob
import time
import json
import sys
import os
import re

# Custom dependencies
import src.dataset as ds
import src.disorder as disorder
import src.msa as msa


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


class Seed(Pipeline):

    # Constructor
    def __init__(
        self, clusters_path, mgnify_path, dask_client,
        mobidb_env=os.environ.copy(), mobidb_cmd=['mobidb_lite.py'],
        muscle_env=os.environ.copy(), muscle_cmd=['muscle']
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
        # Set client (Client(LocalCluster) by default)
        self.dask_client = dask_client

    # Run Seed alignment (for a single cluster)
    def run(self, cluster_names, cluster_dir, verbose=False):

        # # Get clusters output directory
        # cluster_dir = os.path.dirname(cluster_dir)
        # Case output directory does not exits
        if not os.path.exists(cluster_dir):
            # Make directory
            os.mkdir(cluster_dir)

        # Retrieve cluster members from clusters .tsv dataset
        cluster_members = self.get_cluster_members(
            cluster_names=cluster_names,
            verbose=verbose
        )

        # Check cluster members raises error
        self.check_cluster_members(cluster_members)

        # Define sequences accession numbers list
        sequences_acc = [
            # Get i-th sequence accession in cluster named n
            cluster_members[n][i]
            # For each cluster in cluster members
            for n in cluster_members
            # For each component of that cluster
            for i in range(len(cluster_members[n]))
        ]
        # # Debug
        # print(sequences_acc)

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

        # TODO Remove clusters with compositional bias too high

        # Loop through each cluster name to make clusters directories
        for cluster_name in cluster_members:
            # Define cluster subdir
            cluster_subdir = os.path.join(cluster_dir, cluster_name)
            # Check if a subdirectory for current cluster exists
            if not os.path.exists(cluster_subdir):
                # Make sub directory
                os.mkdir(cluster_subdir)

        # Make seed alignment
        self.make_sequences_aln(
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            cluster_dir=cluster_dir
        )

        # # Initialize cluster sequences dict(name: dict(seq acc: fasta entry))
        # cluster_sequences = dict()
        # # Fill cluster sequences dict
        # for cluster_name in cluster_members:
        #     # Add current cluster entry
        #     cluster_sequences.setdefault(cluster_name, dict())
        #     # Loop through each accession number in cluster members
        #     for acc in cluster_members[cluster_name]:
        #         # Add current fasta sequence to current cluster entry
        #         cluster_sequences[cluster_name][acc] = fasta_sequences[acc]

        # # Verbose log
        # if verbose:
        #     print('Retrieved sequences:')
        #     print(json.dumps(cluster_sequences, indent=2))

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

    # Compute compositional bias for every sequence
    def compute_comp_bias(self, cluster_members, fasta_sequences, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None

        # Initialize disorder dict(sequence acc: list of boolean lists)
        disorder_sequences = dict()
        # Get reference to MobiDB Lite instance
        mobidb = self.mobidb
        # Verbose log
        if verbose:
            print('Computing disordered regions for {:d} sequences...'.format(
                len(fasta_sequences)
            ))
        # Run predictior on each  entry in fasta sequence dict(acc: fasta enry)
        futures = self.dask_client.map(
            # Function make a prediction over a single sequence
            lambda acc: mobidb.run(
                sequences={acc: fasta_sequences[acc]},
                verbose=False
            ),
            # Map function over each accession number in fasta sequences
            [*fasta_sequences.keys()]
        )
        # Get distributed computing results
        results = self.dask_client.gather(futures)
        # Loop through every result
        for i in range(len(results)):
            # Get i-th result
            result = results[i]
            # Loop through each sequence accession in i-th result
            for acc in result.keys():
                # Eventually make new entry in disorder dictionart
                disorder_sequences.setdefault(acc, [])
                # Append disorder predictions
                disorder_sequences[acc] += result[acc]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Get execution time summary
            print('Took {:.0f} seconds to compute disordered regions for {:d} fasta sequences'.format(
                time_end - time_beg,  # Time required to distributed search
                len(disorder_sequences)  # Number of fasta sequences
            ))

        # Initialize compossitional bias dict(cluster name: comp bias)
        comp_bias = dict()
        # Verbose log
        if verbose:
            print('Computing compositional bias for {:d} sequences...'.format(
                len(disorder_sequences)
            ))
        # For each cluster, compute compute compositional bias
        for cluster_name in tqdm(cluster_members, disable=(not verbose)):
            # Initialize current cluster disorder predictions
            cluster_disorder = dict()
            # Get disordered regions for current cluster sequences
            for acc in cluster_members[cluster_name]:
                # Store disordered regions for current cluster sequences
                cluster_disorder[acc] = disorder_sequences[acc]
            # Apply a threshold over disorder predictions
            cluster_disorder = disorder.compute_threshold(
                sequences=cluster_disorder,
                threshold=1,
                inclusive=True
            )
            # Compute compositional bias
            comp_bias[cluster_name] = disorder.compute_comp_bias(
                sequences=cluster_disorder
            )

        # Return dict(cluster name: compositional bias)
        return comp_bias

    # Make multiple sequence alignment
    def make_sequences_aln(self, cluster_members, fasta_sequences, cluster_dir, verbose=False):

        # # Initialize dict(cluster name: dict(sequence acc: fasta entry))
        # cluster_sequences = dict()
        # # Fill cluster sequences dict
        # for cluster_name in cluster_members:
        #     # Add current cluster entry
        #     cluster_sequences.setdefault(cluster_name, list())
        #     # Loop through each accession number in cluster members
        #     for acc in cluster_members[cluster_name]:
        #         # Add current fasta sequence to current cluster entry
        #         cluster_sequences[cluster_name] += [fasta_sequences[acc]]

        # Initialize timers
        time_beg, time_end = time.time(), None
        # Verbose log
        if verbose:
            print('Making multiple sequence alignments for {:d} clusters...'.format(
                len(cluster_members)
            ))

        # Get reference to muslce instance
        muscle = self.muscle
        # Feed cluster sequences to alignment algorithm
        futures = self.dask_client.map(
            # Given cluster name, align its inner sequences
            lambda cluster_name: muscle.run(
                sequences=[
                    # Get list of all sequences in current cluster
                    fasta_sequences[acc]
                    # Loop through all current cluster members
                    for acc in cluster_members[cluster_name]
                ],
                acc_regex='>(\S+)',
                verbose=False
            ).to_aln(
                # Define where to store current multiple sequence alignment
                out_path=os.path.join(cluster_dir, cluster_name, 'SEED')
            ),
            # Map function over each cluster name
            [*cluster_members.keys()]
        )
        # Wait for all alignments to finish
        self.dask_client.gather(futures)

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            print('Took {:.0f} seconds to make {:d} multiple sequence alignments'.format(
                time_end - time_beg,
                len(cluster_members)
            ))
