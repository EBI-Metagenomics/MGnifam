# Dependencies
from tempfile import NamedTemporaryFile
import time
import gzip
import stat
import sys
import os
import re


# Check if file is compressed by suffix
def is_gzip(in_path):
    # Return true if file name ends with .gz
    return bool(re.search(r'\.gz$', in_path))


# Open file even if it is gzipped
def open_file(in_path, mode='r', gzip_mode='rt'):
    """Open an eventually compressed file

    Args
    in_path (str)       Path to file which must be opened
    mode (str)          Mode used to open given file path if not compressed
    gzip_mode (str)     Mode used to open given fil if it is compressed

    Return
    (file)              Buffer to opened file with given mode
    """
    # Case file is gzipped
    if is_gzip(in_path):
        # Return uncompressed file buffer
        return gzip.open(in_path, gzip_mode)
    # Otherwise, return common file buffer
    return open(in_path, mode)


# Uncompress file
def gunzip(in_path, out_path=None, out_suffix=''):
    # Open file in read mode
    in_file = open_file(in_path, 'r', 'rt')
    # Case no output path is set
    if not out_path:
        # Set output path as temporary file path
        out_path = NamedTemporaryFile(suffix=out_suffix, delete=False).name
    # Open output path in write mode
    with open(out_path, 'w') as out_file:
        # Loop through every input file path
        for in_line in in_file:
            # Write input line to output file
            out_file.write(in_line)
    # Close input file
    in_file.close()
    # Return path to output file
    return out_path


# Benchmarking function
def benchmark(fn, *args, **kwargs):
    # Define timers
    time_beg, time_end = time.time(), None
    # Execute function with given arguments
    results = fn(*args, **kwargs)
    # Update timers
    time_end = time.time()
    # Return Function results and time taken
    return results, (time_end - time_beg)
