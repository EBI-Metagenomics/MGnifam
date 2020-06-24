"""
Pipeline for making datasets
Takes either mgnify.fa and clusters.tsv files and allows to chunk them in
smalles pieces which can easily be read and searched easily and faster by
leveraging paraellel computation
"""


# Dependencies
import gzip
import os
import re

# Custom dependencies
from src.sequence import Fasta


class Dataset(object):

    # Constructor
    def __init__(self, clusters_path, mgnify_path):
        # Store path to clusters and mgnify
        self.clusters_path = clusters_path
        self.mgnify_path = mgnify_path

    # Split clusters .tsv dataset in chunks
    def chunk_clusters(self, out_path, chunk_size=1e06):
        """Split clusters entries file into chunks
        Given a .tsv[.gz] file, makes chunks of <chunk_size> and stores them
        into out_path directory named according to chunk index.

        Args
        out_path (str)  String containing the path of a generic chunk file,
                        e.g. `chunk{:d}.tsv.gz`. It must be formattable.
        """
        # Get output directory
        out_dir = os.dirname(out_path)
        # Make output directory (chunks will be stored here)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # Open file for reading
        with open_any(self.clusters_path) as in_file:
            # Initialize current chunk content (list of str)
            chunk_content = list()
            # Initialize current line index
            i = 0
            # Loop through every index, line in input file
            for line in in_file:
                # Save current line
                chunk_content.append(line)
                # Case index reached batch size
                if (i+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = i // chunk_size
                    # Persist chunk to disk
                    write_chunk(out_path, chunk_index, chunk_content)
                    # Reinitialize chunk content
                    chunk_content = list()
                # Increase line counter
                i += 1

            # Persist last chunk, if any
            if chunk_content:
                # Define chunk index
                chunk_index = i // chunk_size
                # Persist chunk to disk
                write_chunk(out_path, chunk_index, chunk_content)

    # Split mgnify .fa dataset in chunks
    def chunk_mgnify(self, out_path, chunk_size=1e06):
        """Split MGnify sequences file into chunks
        Given a .fa[.gz] file, makes chunks of <chunk_size> and stores them
        into out_path directory named according to chunk index.
        Note that <chunk_size> is refferred to the number of entries, not the
        number of lines in output chunk, hence chunk sizes are heterogeneous.

        Args
        out_path (str)  String containing the path of a generic chunk file,
                        e.g. `chunk{:d}.fa.gz`. It must be formattable.
        """
        # Get output directory
        out_dir = os.dirname(out_path)
        # Make output directory (chunks will be stored here)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # Open file for reading
        with open_any(self.mgnify_path) as in_file:
            # Initialize current chunk content (list of str)
            chunk_content = list()
            # Initialize current line index
            i = 0
            # Loop through every index, line in input file
            for entry in Fasta.read(in_file):
                # Save current line
                chunk_content.append(entry)
                # Case index reached batch size
                if (i+1) % chunk_size == 0:
                    # Define chunk index
                    chunk_index = i // chunk_size
                    # Persist chunk to disk
                    write_chunk(out_path, chunk_index, chunk_content)
                    # Reinitialize chunk content
                    chunk_content = list()
                # Increase line counter
                i += 1

            # Persist last chunk, if any
            if chunk_content:
                # Define chunk index
                chunk_index = i // chunk_size
                # Persist chunk to disk
                write_chunk(out_path, chunk_index, chunk_content)


# Check if file is gzipped
def is_gzip(path):
    return bool(re.search(r'\.gz$', path))


# Open file, even if gzipped
def open_any(path):
    # Case file is compressed
    if is_gzip(path):
        return gzip.open(path, 'rt')
    # Case file is not compressed
    else:
        return open(path, 'r')


# Write chunk
def write_chunk(path, index, content):
    # Define output path
    chunk_path = path.format(index)
    # Persist data as compressed file
    with gzip.open(chunk_path, 'wt') as chunk_file:
        # Write content
        chunk_file.write('\n'.join(content))
