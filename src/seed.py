"""
Seed alignment pipeline
This script substitutes mgseed and create_seed_from_cluster.pl
"""

# Dependencies
from tqdm import tqdm
import numpy as np
import subprocess
import tempfile
import shutil
import time
import gzip
import re
import os

# Custom dependencies
from src.sequence import Fasta


class Seed(object):

    # Constructor
    def __init__(self, clusters_path, mgnify_path, verbose=False, env=os.environ.copy()):
        # Path to clusters .tsv[.gz] file
        self.clusters_path = clusters_path
        # Path to sequences .fa[.gz] file
        self.mgnify_path = mgnify_path
        # Environmental variables used in external scritps execution
        self.env = env
        # Set verbose output
        self.verbose = verbose

    # Execute pipeline
    def main(self, cluster_name, out_path):

        # Get cluster members (list of squence accession numbers)
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
        # Loop throuth each cluster (sequence accession, fasta) pair
        for acc, fasta in cluster_sequences.items():
            # Case fasta entry is none
            if fasta is None:
                del cluster_sequences[acc]  # Delete entry
        # Case number of retrieved sequences does not match cluster members
        if len(cluster_members) != len(cluster_sequences):
            # Raise new exception
            raise KeyError('Not enough sequences for cluster {:s} found in {:s}'.format(
                cluster_name,  # Cluster name
                self.mgnify_path  # Fasta sequences path
            ))

        # Compute compositional bias (MobiDB Lite)
        comp_bias = self.compute_comp_bias(cluster_sequences)
        # TODO Discard cluster with many low complexity/disorder residues

        # TODO Compute multiple sequence alignment


    # Get cluster members
    def get_cluster_members(self, cluster_name):
        """Retrieve sequences accession numbers
        Takes cluster name as input and searches in given clusters file the
        sequences associated with that file

        Args
        cluster_name (str)  Name of the cluster whose associate sequences
                            accession numbers must be searched

        Return
        (list)              List of sequences accession numbers found to be
                            associated with given cluster name
        """

        # Check if clusters file is compressed
        is_gzip = bool(re.search(r'\.gz$', self.clusters_path))
        # Define an uncompressed file handler
        file_handler = lambda path: open(path, 'r')
        # Case clusters file is compressed
        if is_gzip:
            # Define a compressed file handler
            file_handler = lambda path: gzip.open(path, 'rt')

        # Initialize a list of sequence accession
        sequences_acc = list()
        # Verbose out
        if self.verbose:
            print('Reading clusters file...')
        # Open file with defined file handler
        with file_handler(self.clusters_path) as clusters_file:
            # Loop through every line in file
            for line in tqdm(clusters_file, disable=(not self.verbose)):
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
                sequences_acc.append(found_sequence_acc)
        # Return list of founs sequences accession numbers
        return sequences_acc

    # Get sequences fasta
    def get_sequences_fasta(self, sequences_acc):
        """Retrieve sequences residues
        Takes a list of sequences accessions and search for the associated
        entry by scanning given fasta file headers.

        Args
        sequences_acc (list)    List of sequences accession numbers whose
                                residues must be found in given fasta file

        Return
        (dict(str: str))        Dictionary containing sequences accession
                                numbers as keys and fasta entries as values
        """

        # Check if sequences fasta file is compressed
        is_gzip = bool(re.search(r'\.gz$', self.clusters_path))
        # Define an uncompressed file handler
        file_handler = lambda path: open(path, 'r')
        # Case sequences file is compressed
        if is_gzip:
            # Define a compressed file handler
            file_handler = lambda path: gzip.open(path, 'rt')

        # Initialize output dictionary with sequences accession numbers as keys
        sequences = {acc: None for acc in sequences_acc}
        # Initialize current sequence accession
        curr_acc = None
        # Verbose out
        if self.verbose:
            print('Reading sequences file...')
        # Open file with defined file handler
        with file_handler(self.clusters_path) as clusters_file:
            # Loop through every line in file
            for line in tqdm(clusters_file, disable=(not self.verbose)):
                # Check if line matches fasta header
                match_header = re.search(r'^>([a-zA-Z0-9]+).*[\n\r]*$', line)
                # Case current line is not header
                if not match_header:
                    # Check if there is a valid accession number set
                    if curr_acc is not None:
                        # Store residues in current line
                        sequences[curr_acc] += re.sub(r'[^a-zA-Z]+', '', line)
                    # Otherwise, skip iteration
                # Case current line is header
                else:
                    # Set sequence accession
                    curr_acc = match_header.group(1)
                    # Case current accession is not one of the searched one
                    if curr_acc not in set(sequences.keys()):
                        # Discard sequence accession
                        curr_acc = None
                    # Case current accession must be stored
                    else:
                        # Retrieve header line without special characters
                        header = re.sub(r'[\n\r]*$', '', line)
                        # Store header line (followed by newline)
                        sequences[curr_acc] = header + '\n'
        # Return filled sequences dictionary
        return sequences

    # Get compositional bias for cluster members
    def compute_comp_bias(self, sequences):
        """Compute compositional bias
        Given some sequences fasta entries, compute compositional bias over
        them by running MobiDB Lite classifier.abs($0)

        Args
        sequences (dict(str: str))  Dictionary containing sequence identifier
                                    (e.g. sequence accession number) as key
                                    and associated fasta entry as value

        Return
        (float)                     Compositional bias for given sequences,
                                    i.e. the number of disordered residues
                                    over the total number of residues
        """

        # Define temporary input file
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        # Define temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Make input fasta file
        with open(fasta_file.name, 'w') as ff:
            # Fill FASTA input file
            for sequence_acc in sequences.keys():
                # Write fasta entry
                ff.write(sequences[sequence_acc])

        # Initialize residues dict(sequence acc: np.array(residues))
        residues = {}
        # Fill rediues dictionary
        for sequence_acc in sequences.keys():
            # Split fasta entry in header and residues
            h, r = tuple(sequences[sequence_acc].split('\n'))
            # Turn residues string into numpy array
            r = np.array(list(r))
            # Store residues for current accession
            residues[sequence_acc] = r

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
        # Delete input file
        os.remove(fasta_file)

        # Initialize disorder dict(sequence acc: list((disorder region)))
        # I.e. value is a list of tuples(disordered region begin, end)
        disorder = {acc: list() for acc in sequences.keys()}
        # Get output of MobiDB Lite
        with open(out_file.name, 'r') as of:
            # Loop through every file line
            for line in of:
                # Match output file line format
                match = re.search(r'^(\S+)\s+(\d+)\s+(\d+)', line)
                # Case the format does not match
                if not match:
                    continue  # Skip iteration
                # Get sequence accession
                sequence_acc = match.group(1)
                # Get starting/ending residues of disordered region
                disorder_beg = int(match.group(2))
                disorder_end = int(match.group(3))
                # Initialize sequence accession entry
                disorder.setdefault(sequence_acc, list())
                # Store disordered region begin and end residues
                disorder[sequence_acc] = (disorder_beg, disorder_end)
        # Delete output file
        os.remove(out_file)

        # Initialize number of residues and number of disordered ones
        num_residues = 0
        num_disorder = 0
        # Cast (region begin, end) to np array of shape (#regions, #residues)
        for sequence_acc in disorder.keys():
            # Get number of disordered regions for current accession (rows)
            n = len(disorder[sequence_acc])
            # Get number of residues for current accession (columns)
            m = len(residues[sequence_acc])
            # Generate a zeroes matrix with given shape
            x = np.zeros((n, m), dtype=np.uint8)
            # Set cells in disordered regions to 1
            for i, (beg, end) in enumerate(disorder[sequence_acc]):
                # Note that MobiDB Lite regions begin from 1 and are inclusive
                # Instead, numpy starts from 0 and is exclusive wrt end
                x[i, beg-1:end] = 1
            # Residue is disordered (1) if predicted so at least once
            x = (np.sum(x, axis=0) > 0).astype(np.uint8)
            # Update number of total residues
            num_residues += m
            # Update number of disordered residues
            num_disorder += np.sum(x)

        # Compute compositional bias
        comp_bias = num_disorder / num_residues
        # Return compositional bias
        return comp_bias

    # Compute multiple sequence alignment
    def make_msa(self, sequences, out_path):

        # Initialize temporary input file
        in_file = tempfile.NamedTemporaryFile(delete=False)
        # Initialize temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Fill input fasta file from sequences dict(acc: fasta entry)
        with open(in_file.name, 'w') as in_buffer:
            # Loop through each sequence in input dictionary
            for sequence_acc in sequences.keys():
                # Write fasta entry
                in_buffer.write(sequences[sequence_acc] + '\n')

        # Run Multiple Sequence Alignment script
        _ = subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script actually worked
            encoding='utf-8',  # Output encoding
            env=self.env,  # Environmental variables
            args=[
                'muscle', '-quiet',
                '-in', in_file.name,
                '-out', out_file.name
            ]
        )

        # Persist temporary file to disk
        shutil.copy(out_file.name, out_path + '/' + cluster_name)

        # Delete temporary inpput fasta file
        os.remove(in_file)
        # Delete temporary output fasta file
        os.remove(out_file)


# Unit test
if __name__ == '__main__':

    # Get root
    root = os.path.dirname(os.path.realpath(__file__) + '/..')

    # Define path to clusters file
    CLUSTERS_PATH = '/nfs/production/xfam/pfam/jaina/MGnify_clusters/2019_05/clusters/mgy_seqs.cluster.tsv.gz'
    # Define path to MGnify sequences
    MGNIFY_PATH = '/nfs/production/xfam/pfam/data/mgnify/mgnify.fa.gz'

    # Initialize new seed creation pipeline instance
    seed = Seed(
        clusters_path=CLUSTERS_PATH,
        mgnify_path=MGNIFY_PATH,
        verbose=True
    )

    # Define input cluster name
    cluster_name = 'MGYP000848664103'

    # TEST cluster member retrieval
    # Set beginning time of the method
    time_beg = time.time()
    # Retrieve cluster members
    cluster_members = seed.get_cluster_members(cluster_name)
    # Set ending time of the method
    time_end = time.time()
    # Compute method duration
    time_took = time_end - time_beg
    # Log
    print('Took {:.02f} seconds to retrieve {:d} members of cluster {:s}:'.format(
        time_took,  # Method duration
        len(cluster_members),  # Number of cluster sequences
        cluster_name  # Name of the cluster
    ))
    print(', '.join(cluster_members))  # Plot clusters members)

    # TEST cluster sequences retrieval
    # Reset beginning time of the method
    time_beg = time.time()
    # Retrieve cluster members
    cluster_sequences = seed.get_sequences_fasta(cluster_members)
    # Set ending time of the method
    time_end = time.time()
    # Compute method duration
    time_took = time_end - time_beg
    # Log
    print('Took {:.02f} seconds to retrieve {:d} sequences of cluster {:s}:'.format(
        time_took,  # Method duration
        len(cluster_sequences),  # Number of cluster sequences
        cluster_name  # Name of the cluster
    ))
    print('\n'.join([
        fasta for acc, fasta in cluster_sequences.items()
    ]))

    # Define output fasta file path
    with open(root + '/tmp/test.fasta', 'w') as out_file:
        # Store sequences to file
        for acc, fasta in cluster_sequences.items():
            # Case no entry is empty
            if fasta is None:
                continue  # Skip iteration
            # Otherwise, save line
            out_file.write(fasta + '\n')

    # # Empty cluster sequences
    # # Read cluster sequences from fasta file
    # with open(root + '/tmp/test.fasta', 'r') as in_file:
    #     # Go through every entry in file
    #     for entry in Fasta.read(in_file):
    #         # Split entry in header and residues
    #         header, residues = tuple(entry.split('\n'))
    #         # Get
    #
    # # TEST multiple sequence alignment through MUSCLE
    # # Reset beginning time of the method
    # time_beg = time.time()
    # # Retrieve cluster members
    # cluster_sequences = seed.make_msa(
    #     sequences=cluster_sequences,
    #     out_path=root + '/tmp'
    # )
    # # Set ending time of the method
    # time_end = time.time()
    # # Compute method duration
    # time_took = time_end - time_beg
    # # Log
    # print('Took {:.02f} seconds to make MSA of {:d} sequences of cluster {:s}:'.format(
    #     time_took,  # Method duration
    #     len(cluster_sequences),  # Number of cluster sequences
    #     cluster_name  # Name of the cluster
    # ))
