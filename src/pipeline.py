# Dependencies
from dask.distributed import Client
from glob import glob
import time
import json
import sys
import os
import re

# Custom dependencies
import src.dataset as ds


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
    def __init__(self, clusters_path, mgnify_path, dask_client):
        # Retrieve clusters dataset
        clusters_path = self.check_path(clusters_path)
        self.ds_clusters = ds.Cluster.from_list(clusters_path)
        # Retrieve MGnify dataset
        mgnify_path = self.check_path(mgnify_path)
        self.ds_mgnify = ds.MGnify.from_list(mgnify_path)
        # Set client (Client(LocalCluster) by default)
        self.dask_client = dask_client

    # Run Seed alignment (for a single cluster)
    def run(self, cluster_names):

        # Intialize timers
        time_beg = time.time()
        time_end = None
        time_took = None

        # Debug
        print('Retrieving {:d} clusters:\n{}'.format(
            len(cluster_names),  # Number of clusters
            ', '.join(cluster_names)  # Verbose list of cluster names
        ))
        # Initialize dict(cluster name: list of sequence accession numbers)
        cluster_members = dict()
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
        time_took = time_end - time_beg
        # Debug
        print('Took {:.0f} seconds to search members of {:d} clusters'.format(
            time_took,  # Time required to distributed search
            len(cluster_members.keys())  # Number of clusters
        ))

        # # Debug
        # print(json.dumps(cluster_members))

        # Reset timer
        time_beg = time.time()

        # Get list of clusters with no members
        empty_clusters = [n for n, m in cluster_members.items() if not len(m)]
        # Check if at least one cluster has no member
        if empty_clusters:
            # Define error text
            error_text = 'No sequences found for {:d} clusters: {:s}'.format(
                len(empty_clusters),  # Number of clusters with no members
                ', '.join(empty_clusters)  # Verbose string of empty clusters
            )
            # Raise new error
            raise KeyError(error_text)
        # Free unused memory
        del empty_clusters

        # Initialize dict(sequence accession: fasta entry)
        cluster_sequences = dict()
        # Submit jobs for searching cluster members sequences
        futures = self.dask_client.map(
            # Search for cluster members in current chunk
            lambda ds: ds.search([
                # Get i-th sequence accession in cluster named n
                cluster_members[n][i]
                for n in cluster_members.keys()
                for i in range(len(cluster_members[n]))
            ]),
            # List of cluster datasets
            self.ds_mgnify
        )
        # Retrieve found sequences (list of dict)
        results = self.dask_client.gather(futures)
        # # Debug
        # print('Debug', json.dumps(results))
        # Loop through resulting sequences dictionaries
        for i in range(len(results)):
            # Loop through retrieved sequence accession numbers
            for sequence_acc in results[i]:
                # Update cluster sequences dictionary
                cluster_sequences[sequence_acc] = results[i][sequence_acc]

        # Update timers
        time_end = time.time()
        time_took = time_end - time_beg
        # Debug
        print('Took {:.0f} seconds to search for {:d} sequences:'.format(
            time_took,  # Time required to distributed search
            len(cluster_sequences)  # Number of sequences
        ))
        # print(json.dumps(cluster_sequences))

        # Initialize list of clusters with mismatching sequences
        mismatch_clusters = list()
        # Fill mismatching clusters
        for cluster_name in cluster_members.keys():
            # Get sequence accession numbers for current cluster
            sequences_acc = set(cluster_members[cluster_name])
            # Case some sequence has not been retrieved
            if sequences_acc - set(cluster_sequences.keys()):
                # Add cluster name to list of mismatching ones
                mismatch_clusters.append(cluster_name)
        # Case there is at least one mismatching cluster
        if mismatch_clusters:
            # Define new error message
            error_text = 'Could not retrieve all sequences for {:d} clusters: {}'.format(
                len(mismatch_clusters),  # Number of mismatching clusters
                ', '.join(mismatch_clusters)  # Mismatching clusters
            )
            # Raise new error
            raise KeyError(error_text)
        # Free unused memory
        del mismatch_clusters

        # Debug
        print('Retrieved sequences:')
        print(json.dumps(cluster_sequences, indent=2))

    # Retrieve cluster members for a batch of clusters
    def get_cluster_members(self, cluster_names):
        raise NotImplementedError

    # Retrieve cluster sequences for given sequences accession numbers
    def get_fasta_sequences(self, sequences_acc):
        raise NotImplementedError

    # Compute compositional bias for every sequence
    def compute_comp_bias(self, fasta_sequences):
        raise NotImplementedError

    # Make multiple sequence alignment
    def align_sequences(self, fasta_sequences):
        raise NotImplementedError
