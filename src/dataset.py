# Dependencies
import glob
import gzip
import sys
import os
import re

# Custom dependencies
from src.sequence import Fasta


class Dataset(object):

    # Constructor
    def __init__(self, path):
        # Store path to dataset
        self.path = path

    # Abstract
    def to_chunk(self, *args, **kwargs):
        raise NotImplementedError

    # Abstract
    def search(self, *args, **kwargs):
        raise NotImplementedError

    @classmethod
    def from_list(cls, paths):
        return [cls(path) for path in paths]

    @classmethod
    def from_str(cls, path):
        return cls.from_list(glob.glob(path))

    @staticmethod
    def is_gzip(path):
        # Return true if file name ends with .gz
        return bool(re.search(r'\.gz$', path))

    @staticmethod
    def open_file(path):
        # Case file is compressed
        return gzip.open(path, 'rt') if Dataset.is_gzip(path) else open(path, 'r')

    @staticmethod
    def write_chunk(path, index, content, sep=''):
        # Define output path
        chunk_path = path.format(index)
        # Persist data as compressed file
        with gzip.open(chunk_path, 'wt') as chunk_file:
            # Write content
            chunk_file.write(sep.join(content))


class Mgnify(Dataset):

    # Chunking function
    def to_chunks(self, chunk_path='chunk{:03d}.tsv.gz', chunk_size=1e07):
        """Split MGnify sequences file into chunks
        Given a .fa[.gz] file, makes chunks of <chunk_size> and stores them
        into out_path directory named according to chunk index.
        Note that <chunk_size> is refferred to the number of entries, not the
        number of lines in output chunk, hence chunk sizes are heterogeneous.

        Args
        chunk_path (str)    String containing the path of a generic chunk file,
                            e.g. `chunk{:d}.fa.gz` (must be formattable)
        chunk_size (int)    Maximum number of fasta entries to be stored in
                            each chunk
        """
        # Get output directory
        chunks_dir = os.path.dirname(chunk_path)
        # Case given output directory does not exist
        if not os.path.exists(chunks_dir):
            # Attempt to make a new output directory
            os.mkdir(chunks_dir)
        # Initialize current chunk (batch of fasta sequences entries)
        seq_batch = list()
        # Initialize sequence index
        seq_index = 0
        # Open file for reading
        with self.open_file(self.path) as file:
            # Loop through every index, line in input file
            for entry in Fasta.read(file):
                # Save current line
                seq_batch.append(entry)
                # Case index reached batch size
                if (seq_index+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = seq_index // chunk_size
                    # Persist chunk to disk
                    self.write_chunk(chunk_path, chunk_index, seq_batch, sep='\n')
                    # Reinitialize chunk content
                    seq_batch = list()
                # Increase line counter
                seq_index += 1
            # Persist last chunk, if any
            if seq_batch:
                # Define chunk index
                chunk_index = seq_index // chunk_size
                # Persist chunk to disk
                self.write_chunk(chunk_path, chunk_index, seq_batch, sep='\n')

    # Search function
    def search(self, sequences_acc):
        raise NotImplementedError


class Cluster(Dataset):

    # Chunking function
    def to_chunks(self, chunk_path='chunk{:03d}.fa.gz', chunk_size=1e07):
        """Split clusters entries file into chunks
        Given a .tsv[.gz] file, makes chunks of <chunk_size> and stores them
        into out_path directory named according to chunk index.

        Args
        chunk_path (str)    String containing the path of a generic chunk file,
                            e.g. `chunk{:d}.tsv.gz` (must be formattable)
        chunk_size (int)    Maximum number of lines to be stored in each chunk
        """
        # Get output directory
        chunks_dir = os.path.dirname(chunk_path)
        # Case given output directory does not exist
        if not os.path.exists(chunks_dir):
            # Attempt to make a new output directory
            os.mkdir(chunks_dir)
        # Initialize current chunk (batch of sequence accession numbers)
        seq_batch = list()
        # Initialize sequence index
        seq_index = 0
        # Open input file
        with self.open_file(self.path) as file:
            # Loop through every index, line in input file
            for line in file:
                # Save current line
                seq_batch.append(line)
                # Case index reached batch size
                if (seq_index+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = seq_index // chunk_size
                    # Persist chunk to disk
                    self.write_chunk(chunk_path, chunk_index, seq_batch)
                    # Reinitialize chunk content
                    seq_batch = list()
                # Increase line counter
                seq_index += 1
            # Persist last chunk, if any
            if seq_batch:
                # Define chunk index
                chunk_index = seq_index // chunk_size
                # Persist chunk to disk
                self.write_chunk(chunk_path, chunk_index, seq_batch)

    # Search function
    def search(self, cluster_names, verbose=False):
        """Retrieve sequences accession numbers
        Takes list of cluster names as input and searches in given clusters
        file the sequences accession numbers associated with that cluster name

        Args
        cluster_names (list(str))   Names of the clusters whose sequences
                                    members must be retrieved
        verbose (bool)              Wether to print out verbose log

        Return
        (dict(str:list))            Dictionary whose keys are cluster names and
                                    values are lists of strings containing
                                    sequences accession numbers associated with
                                    the key's cluster name
        """
        # Parse cluster names to set
        cluster_names = set(cluster_names)
        # Initialize output dictionary
        sequences_acc = dict(cluster: list() for cluster in cluster_names)
        # Open dataset file
        with self.open_file(self.path) as file:
            # Define iterator
            line_iterator = tqdm(
                clusters_file,  # File line iterator
                disable=(not verbose),  # Set verbose
                file=sys.stdout  # Force printing to stdout
            )
            # Loop through every line in file
            for line in line_iterator:
                # Match cluster name and sequence accession
                match = re.search(r'^([a-zA-Z0-9]+)[ \t]+([a-zA-Z0-9]+)', line)
                # Case the line format does not match
                if not match:
                    continue  # Skip iteration
                # Retrieve cluster name and sequence accession
                found_cluster_name = str(match.group(1))
                found_sequence_acc = str(match.group(2))
                # Case cluster name does not match the searched one
                if found_cluster_name not in cluster_names:
                    continue  # Skip iteration
                # Otherwise, store found sequence accesssion
                sequences_acc[found_cluster_name].append(found_sequence_acc)
        # Return dict (cluster name: sequences accession numbers)
        return sequences_acc
