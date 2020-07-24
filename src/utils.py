# Dependencies
from tempfile import NamedTemporaryFile
from glob import glob
import time
import gzip
import stat
import sys
import os
import re


def get_paths(in_paths):
    # Check input path type: string
    if isinstance(in_paths, str):
        # Use glob to find files referenced by unix string
        return glob(in_paths)
    # Case input path is already a list
    elif isinstance(in_paths, list):
        # Return given list
        return in_paths
    # Case input type is not valid
    raise ValueError('Given path is not valid')


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
    return (time_end - time_beg), results


# Turn memory string into bits
def as_bytes(memory_str):
    # Memory multipliers
    multipliers = {'': 1, 'B': 1, 'KB': 1e03, 'MB': 1e06, 'GB': 1e09, 'TB': 10e12}
    # Remove trailing spaces
    memory_str = memory_str.strip()
    # Check that input string matches format
    match = re.search(r'^(\d+\.{,1}|\d+\.\d+|\.\d+)\s*(\S*)$', memory_str)
    # # Debug
    # print('Matches?', bool(match))
    # Case string matches
    if match:
        # Otherwise, get either numeric value and multiplier
        value = float(match.group(1))
        size = str(match.group(2)).upper()
        # # Debug
        # print('Value:', value)
        # print('Size:', size)
        # Case size mathces a multiplier
        if size in set(multipliers.keys()):
            # Return float value times multiplier
            return value * multipliers[size]
    # Raise new error
    raise ValueError('Wrong memory string input format')
