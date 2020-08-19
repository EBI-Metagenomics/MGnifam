# Dependencies
from src.sequences import fasta_iter
from src.utils import open_file
from tqdm import tqdm
import glob
import sys
import os
import re


class Dataset(object):

    # Constructor
    def __init__(self, path):
        # Store path to dataset
        self.path = path

    # Just a wrapper for get length
    def __len__(self):
        # Call get length method
        return self.get_length()

    # Retrieve dataset length
    def get_length(self):
        raise NotImplementedError

    # (Abstract) Writes a dataset partition to file
    def to_chunk(self, *args, chunk_path='chunk{:d}', chunk_size=1e06, **kwargs):
        raise NotImplementedError

    # (Abstract) Retrieves entries according to some parameters
    def search(self, *args, ret_length=False, verbose=False, **kwargs):
        raise NotImplementedError

    @classmethod
    def from_list(cls, paths):
        return [cls(path) for path in paths]

    @classmethod
    def from_str(cls, path):
        return cls.from_list(glob.glob(path))

    @staticmethod
    def write_chunk(path, index, content, sep=''):
        # Define output path
        chunk_path = path.format(index)
        # Persist data as compressed file
        with open_file(chunk_path, 'w', 'wt') as chunk_file:
            # Write content
            chunk_file.write(sep.join(content))


class Fasta(Dataset):

    # Get length (number of fasta sequences)
    def get_length(self):
        # Initialize output length
        length = 0
        # Open underlying file
        with open_file(self.path) as file:
            # Loop through each entry in input fasta file
            for entry in fasta_iter(file):
                # Update dataset length
                length += 1
        # Return dataset length
        return length

    # Get longest sequence
    def get_longest(self):
        # Initialize current longest entry and its length (number of residues)
        longest_seq, longest_len = '', 0
        # Initailize number of sequences
        num_sequences = 0
        # Open inner dataset file path
        with open_file(self.path) as file:
            # Loop through each file entry
            for entry in fasta_iter(file):
                # Split current entry in header and residues
                header, residues = tuple(entry.split('\n'))
                # Get current sequence and its number of residues
                curr_seq, curr_len = entry, len(residues)
                # Case current sequence is longer than longest
                if curr_len > longest_len:
                    # Update longest sequence and its length
                    longest_seq, longest_len = curr_seq, curr_len
                # Unpate number of sequences
                num_sequences += 1
        # Return either longest sequence and its length
        return longest_seq, longest_len, num_sequences

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

        Raise
        (FileNotFoundError) If given chunk path is not valid
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
        with open_file(self.path) as file:
            # Loop through every index, line in input file
            for entry in fasta_iter(file):
                # Save current line
                seq_batch.append(entry)
                # Case index reached batch size
                if (seq_index+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = int(seq_index // chunk_size)
                    # Persist chunk to disk
                    self.write_chunk(chunk_path, chunk_index, seq_batch, sep='\n')
                    # Reinitialize chunk content
                    seq_batch = list()
                # Increase line counter
                seq_index += 1
            # Persist last chunk, if any
            if seq_batch:
                # Define chunk index
                chunk_index = int(seq_index // chunk_size)
                # Persist chunk to disk
                self.write_chunk(chunk_path, chunk_index, seq_batch, sep='\n')
        # Define number of chunks
        num_chunks = chunk_index + 1
        # Return number of chunks
        return num_chunks

    # Search function
    def search(self, sequences_acc, ret_length=False, verbose=False):
        """Retrieve sequences residues
        Takes a list of sequences accessions and search for the associated
        entry by scanning underlying fasta file headers.

        Args
        sequences_acc (list)    List of sequences accession numbers whose
                                residues must be found in given fasta file
        ret_length (bool)       Wether to return the length of the searched
                                target dataset (disables early stopping
                                criterion)
        verbose (bool)          Whether to print out verbose log

        Return
        (dict(str: str))        Dictionary containing sequences accession
                                numbers as keys and fasta entries as values
        """
        # Cast cluster names to set
        sequences_acc = set(sequences_acc)
        # Initialize output dict(sequence acc: fasta entry) and length
        sequences, length = dict(), 0
        # Verbose out
        if verbose:
            print('Reading sequences file', self.path)
        # Open file with defined file handler
        with open_file(self.path) as file:
            # Define fasta entries iterator
            tqdm_iter = tqdm(
                fasta_iter(file),  # Input iterator
                disable=(not verbose),  # Set verbose
                file=sys.stdout  # Force printing to stdout
            )
            # Loop through each entry in input fasta file
            for entry in tqdm_iter:
                # Split entry in header and residues
                header, resiudes = entry.split('\n')
                # Get accession number from header
                acc = re.search(r'^>(\S+)', header).group(1)
                # Case accession is one of the searched ones
                if acc in sequences_acc:
                    # Store entry
                    sequences[acc] = entry
                # Case all sequences have been found
                if (not ret_length) and (len(sequences) == len(sequences_acc)):
                    break  # Early stopping
                # Case length must be returned
                elif ret_length:
                    length += 1
        # Case length must be returned
        if ret_length:
            return sequences, length
        # Case only sequences must be returned
        else:
            return sequences


class LinClust(Dataset):

    # Chunking function
    def to_chunks(self, chunk_path='chunk{:03d}.fa.gz', chunk_size=1e07):
        """Split clusters entries file into chunks
        Given a .tsv[.gz] file, makes chunks of <chunk_size> and stores them
        into out_path directory named according to chunk index.

        Args
        chunk_path (str)    String containing the path of a generic chunk file,
                            e.g. `chunk{:d}.tsv.gz` (must be formattable)
        chunk_size (int)    Maximum number of lines to be stored in each chunk

        Raise
        (FileNotFoundError) If given chunk path is not valid
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
        with open_file(self.path) as file:
            # Loop through every index, line in input file
            for line in file:
                # Save current line
                seq_batch.append(line)
                # Case index reached batch size
                if (seq_index+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = int(seq_index // chunk_size)
                    # Persist chunk to disk
                    self.write_chunk(chunk_path, chunk_index, seq_batch)
                    # Reinitialize chunk content
                    seq_batch = list()
                # Increase line counter
                seq_index += 1
            # Persist last chunk, if any
            if seq_batch:
                # Define chunk index
                chunk_index = int(seq_index // chunk_size)
                # Persist chunk to disk
                self.write_chunk(chunk_path, chunk_index, seq_batch)
        # Define number of chunks
        num_chunks = chunk_index + 1
        # Return number of chunks
        return num_chunks

    # Search function
    def search(self, cluster_names, verbose=False):
        """Retrieve sequences accession numbers
        Takes list of cluster names as input and searches in given clusters
        file the sequences accession numbers associated with that cluster name

        Args
        cluster_names (list(str))   Names of the clusters whose sequences
                                    members must be retrieved
        verbose (bool)              Whether to print out verbose log

        Return
        (dict(str:list))            Dictionary whose keys are cluster names and
                                    values are lists of strings containing
                                    sequences accession numbers associated with
                                    the key's cluster name
        """
        # Cast cluster names to set
        cluster_names = set(cluster_names)
        # Initialize output dictionary
        sequences_acc = {cluster: list() for cluster in cluster_names}
        # Verbose out
        if verbose:
            print('Reading clusters file', self.path)
        # Open dataset file
        with open_file(self.path) as file:
            # Define iterator
            line_iterator = tqdm(
                file,  # File line iterator
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
