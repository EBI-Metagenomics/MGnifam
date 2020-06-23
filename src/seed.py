"""
Seed alignment pipeline
This script substitutes mgseed and create_seed_from_cluster.pl
"""

# Dependencies
import subprocess
import tempfile
import gzip
import re
import os


class Seed(object):

    # Constructor
    def __init__(self, clusters_path, mgnify_path, out_path, env=os.environ.copy()):
        # Path to clusters .tsv[.gz] file
        self.clusters_path = clusters_path
        # Path to sequences .fa[.gz] file
        self.mgnify_path = mgnify_path
        # Path to output directory
        self.out_path = out_path
        # Environmental variables used in external scritps execution
        self.env = env

    # Wrapper for main
    def __call__(self, *args, **kwargs):
        # Execute main
        return self.main(*args, **kwargs)

    # Execute pipeline
    def main(self, cluster_name):

        # Initialize cluster sequences dict(accession: fasta)
        cluster_sequences = dict()

        # Get cluster members (squence accession numbers)
        cluster_members = self.get_cluster_members(cluster_name)
        # Case there is no cluster member available
        if not cluster_members:
            # Raise new exception
            raise KeyError('{:s} cluster not found in {:s}'.format(
                cluster_name,  # Cluster name
                self.clusters_path  # Clusters file path
            ))

        # Get fasta sequences
        cluster_sequences = self.get_sequences_fasta(cluster_members)
        # Case number of retrieved sequences does not match cluster members
        if len(cluster_members) != len(cluster_sequences):
            # Raise new exception
            raise KeyError('Not enough sequences for cluster {:s} found in {:s}'.format(
                cluster_name,  # Cluster name
                self.mgnify_path  # Fasta sequences path
            ))

        # TODO Compute compositional bias (MobiDB Lite)
        # TODO Discard clusters with high low complexity/disorder regions

        # TODO Compute multiple sequence alignment

    # Get cluster members
    def get_cluster_members(self, cluster_name):

        # Check if clusters file is compressed
        is_gzip = bool(re.search(r'\.gz$', self.clusters_path))
        # Define an uncompressed file handler
        file_handler = lambda path: open(path, 'r')
        # Case clusters file is compressed
        if is_gzip:
            # Define a compressed file handler
            file_handler = lambda path: gzip.open(path, 'rt')

        # Initialize a list of sequence accession
        cluster_members = list()
        # Open file with defined file handler
        with file_handler(self.clusters_path) as clusters_file:
            # Loop through every line in file
            for line in clusters_file:
                # Match cluster name and sequence accession
                match = re.search(r'^([a-zA-Z0-9]+)[ \t]+([a-zA-Z0-9]+)', line)
                # Case the line format does not match
                if not match:
                    continue  # Skip iteration
                # Retrieve cluster name and sequence accession
                found_cluster_name = str(match.group(1))
                found_sequence_acc = str(match.group(2))
                # Case cluster name does not match the searched one
                if found_cluster_name != cluster_name:
                    continue  # Skip iteration
                # Otherwise, store found sequence accesssion
                cluster_members.append(found_sequence_acc)
        # Return list of cluster members
        return cluster_members

    # Get sequences fasta
    def get_sequences_fasta(self, cluster_members):

        # Check if sequences fasta file is compressed
        is_gzip = bool(re.search(r'\.gz$', self.clusters_path))
        # Define an uncompressed file handler
        file_handler = lambda path: open(path, 'r')
        # Case sequences file is compressed
        if is_gzip:
            # Define a compressed file handler
            file_handler = lambda path: gzip.open(path, 'rt')

        # Cast list of cluster members to set
        cluster_members = set(cluster_members)

        # Initialize a list of sequence accession
        cluster_sequences = {'': ''}
        # Open file with defined file handler
        with file_handler(self.clusters_path) as clusters_file:
            # Loop through every line in file
            for line in clusters_file:
                # Check if line matches fasta header
                match_header = re.search(r'^>(.*)[\n\r]*$', line)
                # Define previous header
                curr_header = [*sequences.keys()][-1]
                # Case current line is header
                if match_header:
                    # Define key (remove special characters)
                    header = match_header.group(1)
                    # Get sequence accession
                    curr_acc = re.search(r'^([a-zA-Z0-9]+)[ \t]+', curr_header)
                    # Delete previous header if not cluster member
                    if curr_acc not in cluster_members:
                        # Delete previous header
                        del sequences[curr_header]
                    # Update current header
                    curr_header = header
                    # Define new entry in sequence dictionary
                    sequences[curr_header] = ''
                # Case current line is not header
                else:
                    # Append current sequence to its header
                    sequences[curr_header] += re.sub(r'[^a-zA-Z]+', '', line)
            # Get last sequence header
            curr_header = [*cluster_sequences.keys()][-1]
            # Get last sequence accession
            curr_acc = re.search(r'^([a-zA-Z0-9]+)[ \t]+', curr_header)
            # Delete last sequence if not cluster member
            del cluster_sequences[curr_header]
        # Return retrieved sequences
        return cluster_sequences

    # Get compositional bias for cluster members
    def get_compositional_bias(self, sequences):
        # Define temporary input file
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        # Define temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Add header entry
        with open(fasta_file.name, 'w') as ff:
            # Fill FASTA input file
            for header in sequences.keys():
                # Write header row
                ff.write('>' + header + '\n')
                # Write content row
                ff.write(sequences[header] + '\n')

        # Run MobiDB Lite
        _ = subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script actually worked
            encoding='utf-8',  # Output encoding
            env=self.env,  # Environmental variables
            # Executable arguments
            args=[
                'python3', 'mobidb_lite.py', '-o', out_file, fasta_file.name
            ]
        )

        # Get output of MobiDB Lite
        with open(out_file.name, 'r') as of:
            # Loop through every file line
            for line in of:
                # Match output file line format
                match = re.search(r'^([a-zA-Z0-9]+)[ \t]+(\d+)[ \t]+(\d+)', line)
                # Case the format does not match
                if not match:
                    continue  # Skip iteration
                # Get sequence accession
                sequence_acc = match.group(1)
                # Get starting/ending residues of disordered region
                disorder_beg, disorder_end = int(match.group(2)), int(match.group(3))
                # Compute number of disordered residues
                num_disorder = disorder_end - disorder_beg

        # Delete input file
        os.remove(fasta_file)
        # Delete output file
        os.remove(out_file)
