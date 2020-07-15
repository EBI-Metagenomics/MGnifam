"""HMM scan and HMM build

This library wraps HMMER functions for either building an Hidden Markov Model
(HMM) starting from a Multiple Sequence Alignment (MSA) and for searching a
HMM against a given dataset (e.g. UniProt).
"""


# Dependecies
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
from itertools import product
from tqdm import tqdm
import numpy as np
import subprocess
import traceback
import glob
import time
import math
import json
import sys
import os
import re


# Custom dependencies
import src.dataset as ds
from src.utils import open_file, gunzip, benchmark
from src.msa import MSA


class HMM(object):

    # Constructor
    def __init__(self, name='', length=0, alphabet='amino'):
        # HMM attributes
        self.name = name
        self.length = length
        self.alphabet = alphabet

    @classmethod
    def from_file(cls, path):
        # Define new dict containing default HMM parameters
        params = {'name': '', 'length': 0, 'alphabet': 'amino'}
        # Open input file
        file = open_file(path, 'r')
        # Loop through each line in file
        for line in file:
            # Check if line is HMM name
            is_name = re.search(r'^NAME\s+(\S+)', line)
            # Case line is HMM name
            if is_name:
                # Store name
                params['name'] = str(is_name.group(1))
                # Go to next iteration
                continue
            # Check if line is HMM length
            is_length = re.search(r'^LENG\s+(\d+)', line)
            # Case line is HMM length
            if is_length:
                # Store HMM length
                params['length'] = int(is_length.group(1))
                # Go to next iteration
                continue
            # Check if line is HMM alphabet
            is_alphabet = re.search(r'^ALPH\s+(\S+)', line)
            # Case line is HMM alphabet
            if is_alphabet:
                # Store HMM alphabet
                params['alphabet'] = str(is_alphabet.group(1))
                # Go to next iteration
                continue
        # Close file buffer
        file.close()
        # Return new HMM instance with retrieved params
        return cls(**params)

    @classmethod
    def concat(cls, in_paths, out_path):
        """Create HMM library
        Create an HMM library by concatenating multiple HMM files

        Args
        in_paths (list)     List of HMM files/libraries to concatenate
        out_path (str)      Path to resulting HMM library file

        Raise
        (FileNotFoundError) Error when trying to open one input file or when
                            creating the output file
        """
        # Open output file
        with open(out_path, 'w') as out_file:
            # Loop through each input path
            for in_path in in_paths:
                # Open input path
                in_file = open(in_path, 'r')
                # Loop through every line in input file
                for line in in_file:
                    # Write out line in output file
                    out_file.write(line)
                # Close input file
                in_file.close()

class HMMER(object):

    def __init__(self, cmd, env=os.environ.copy()):
        self.cmd = cmd
        self.env = env

    def run(self, *args, **kwargs):
        raise NotImplementedError

    # Just a wrapper for run(...) method
    def __call__(self, *args, **kwargs):
        return self.run(*args, **kwargs)


class HMMBuild(HMMER):

    # Set available alphabets
    alphabet = set(['amino', 'dna', 'rna'])

    # Constructor
    def __init__(
        self,
        cmd=['/nfs/production/xfam/pfam/software/bin/hmmbuild'],
        env=os.environ.copy()
    ):
        # Call parent constructor
        super(HMMBuild, self).__init__(cmd=cmd, env=env)

    # Run HMM build script
    def run(
        self, msa_path, out_path, name=None, log_path='/dev/null',
        alphabet='amino'
    ):

        # Initialize command line
        cmd = [*self.cmd]

        # Set HMM name (query name in hmmsearch)
        cmd += ['-n', name] if name else []

        # Case alphabet type is not valid
        if alphabet not in self.alphabet:
            # Raise new error
            raise ValueError('Worng value for alphabet: {}'.format(alphabet))
        # Set type of alphabet used
        cmd += ['--{:s}'.format(alphabet)]

        # Set log path (default /dev/null, i.e. none)
        cmd += ['-o', log_path]

        # Set output file path
        cmd += [out_path]
        # Set input MSA file in script
        cmd += [msa_path]

        # # Define new input MSA temporary file
        # aln_file = tempfile.NamedTemporaryFile(suffix='.aln', delete=False)
        # # Write out input MSA
        # msa.to_aln(aln_file.name)

        # Run HMM build, return CompletedProcess instance
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )


class HMMSearch(HMMER):

    # Constructor
    def __init__(
        self,
        cmd=['/nfs/production/xfam/pfam/software/bin/hmmsearch'],
        env=os.environ.copy()
    ):
        # Call parent constructor
        super(HMMSearch, self).__init__(cmd=cmd, env=env)

    def run(
        self, hmm_path, ds_path, out_path='/dev/null', cpu=None, z=None,
        # Other tabular output formats
        tblout_path='', domtblout_path='', pfamtblout_path=''
    ):
        """Run hmmsearch with some options

        Args
        hmm_path (str)          Path to Hidden Markov Model (HMM) previously
                                created through hmmbuild
        ds_path (str)           Path to dataset against which given HMM must be
                                searched
        out_path (str)          Path to output file
        cpu (int)               Number of cpus to be used in multithreading
        z (int)                 Z value (length of the dataset) to be used
                                in E-value computation
        tblout (str)            Path to tabular per-sequence output
        domtblout (str)         Path to tabular per-domain output
        pfamtblout (str)        Path to pfam-formatted tabular output

        Return
        (CompletedProcess)      Instance returned by subprocess.run(...) for a
                                successfully completed process (contains stdout)

        Raise
        (CalledProcessError)    Raised in case hmmbuild process fails
        """
        # Set hmmsearch executable
        cmd = [*self.cmd]

        # Set output path, if any
        cmd += ['-o', out_path] if out_path else []

        # Set number of cpus used, if any
        cmd += ['--cpu', str(cpu)] if cpu else []

        # Set Z value
        cmd += ['-Z', str(z)] if z else []

        # Set tabular formats
        cmd += ['--tblout', tblout_path] if tblout_path else []
        cmd += ['--domtblout', domtblout_path] if domtblout_path else []
        cmd += ['--pfamtblout', pfamtblout_path] if pfamtblout_path else []

        # Set input hmm file
        cmd += [hmm_path]
        # Set input dataset path
        cmd += [ds_path]

        # Run HMM search
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )

    @staticmethod
    def get_resources(max_mem, model_len, max_cpus=1, long_seq=4e04, num_bytes=48):
        """Computes the resources needed to run hmmscan

        Number of cpus is the maximum number of cpus satisfying the following:
        max     num_cpus
        s.t.    (model_len * long_seq * num_bytes * num_cpus) / 1e09 > mem_threshold
        where   mem_threshold is the maximum amount of memory that can be allocated,
                model_len is the length of the model (found in HMM file),
                long_seq is the length of the longest sequence compared with HMM,
                num_bytes is the number of bytes in the dp

        Args
        max_mem (int)           Maximum number of Gb that can be allocated to
                                run hmmscan
        long_seq (int)          Length of the longest sequence in hmmscan target
                                dataset
        model_len (int)         Length of the HMM model
        max_cpus (int)          Maximum number of available CPUs
        num_bytes (int)         Number of bytes in DP

        Return
        (int)                   Number of cpus needed
        (int)                   Number of Gb of memory needed

        Raise
        (MemoryError)           When the memory threshold is exceeded in any case
        """
        # # Compute threshold
        # threshold = (max_mem * 1e09) / (model_len * long_seq * num_bytes)
        # # Debug
        # print('Threshold =', threshold)
        # print('  longest sequence =', int(long_seq))
        # print('  Gb =', int(1e09))
        # # Case threshold is above maximum number of cpus
        # if max_cpus < threshold:
        #     # Raise new memory exception
        #     raise MemoryError('HMMScan exceeds maximum memory threshold')
        # # Otherwise, loop from maximum number of cpus available down to 1
        # for num_cpus in range(max_cpus, 0, -1):
        #     # Case current number of cpus does not satisfy threshold
        #     if num_cpus < threshold:
        #         continue  # Skip iteration
        #     # Compure required memory
        #     req_memory = int((model_len * long_seq * num_bytes * num_cpus) / 1e09)
        #     # Otherwise, return current number of cpus
        #     return num_cpus, req_memory
        # TODO Check if minimum number of CPUs does not exceed memory threshold
        min_mem = math.ceil(model_len * long_seq * num_bytes * 1) / 1e09
        # Case minimum memory exceeds memory threshold
        if min_mem > max_mem:
            # Raise new memory exception
            raise MemoryError('hmmscan exceeds maximum memory threshold')
        # Loop until required memory does not exceed maximum memory threshold
        for num_cpus in range(max_cpus, 0, -1):
            # Compure required memory
            req_mem = int((model_len * long_seq * num_bytes * num_cpus) / 1e09)
            # Debug
            print('Required resources:')
            print('  number of CPUs:', num_cpus)
            print('  required Gb of memory:', req_mem)
            # Case required memoryexceed maximum threshold
            if req_mem > max_mem:
                continue # Go to next iteration
            # Otherwise, return current number of cpus
            return num_cpus, req_mem

# domtblout iterator
def iter_domtblout(iterable):
    """domtblout generator

    Args
    iter (iterable)     Iterable which returns a domtblout row as a string at
                        each iteration

    Return
    (generator)         Generator returning each line of the domtblout table as
                        a list of tuples, whose first value is the column name
                        and the second is the actual value of the cell
    """
    # Skip first line
    _ = next(iterable)
    # Save second line (header, contains column names)
    header_line = re.sub(r'[\n\r]', '', next(iterable))
    header_line = re.sub(r'^#', ' ', header_line)  # Handle first character
    # Save third line (bounds, define wether a column starts and ends)
    bounds_line = re.sub(r'[\n\r]', '', next(iterable))
    bounds_line = re.sub(r'^#', '-', bounds_line)  # Handle first character

    # Initialize columns names and boundaries lists
    columns, bounds = list(), list()
    # Loop through each character in bounds line
    for j in range(len(bounds_line)):
        # Define if current character is gap
        curr_gap = bool(re.search(r'\s', bounds_line[j]))
        # Case this is the first character
        if (j==0) and (not curr_gap):
            # Define start of the first column
            bounds.append((0, None))
        # Otherwise
        else:
            # Define if previous character is gap
            prev_gap = bool(re.search(r'\s', bounds_line[j-1]))
            # Case current character is gap while previous is not
            if curr_gap and (not prev_gap):
                # Update last column boundary
                bounds[-1] = (bounds[-1][0], j)
            # Case current character is not gap, while previous is
            elif (not curr_gap) and (prev_gap):
                # Add new column entry
                bounds.append((j, None))

    # Update last bound to go until the end of the line
    bounds[-1] = (bounds[-1][0], None)

    # Define column names
    columns = [(header_line[b[0]:b[1]]).strip() for b in bounds]

    # Go through each remaining line in file
    for line in iterable:
        # Remove newlines
        line = re.sub(r'[\n\r]', '', line)
        # Check if line is comment
        is_comment = re.search(r'^#', line)
        # Case current line is comment
        if is_comment:
            continue  # Skip iteration
        # Yield cureent output row
        yield [
            # Tuple (column name, cell value)
            (columns[j], (line[bounds[j][0]:bounds[j][1]]).strip())
            # Loop through each column
            for j in range(len(bounds))
        ]

# TODO tblout iterator
def iter_tblout(iterable):
    raise NotImplementedError

# TODO pfamtblout iterator
def iter_pfamtblout(iterable):
    raise NotImplementedError


# Unit testing
if __name__ == '__main__':

    # Test-only dependencies
    from dask_jobqueue import LSFCluster
    from dask.distributed import Client

    # Define project root path
    ROOT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
    # Define path to temporary folder
    TMP_PATH = os.path.join(ROOT_PATH, 'tmp')
    # Path to test seed file
    SEED_PATH = os.path.join(TMP_PATH, 'MGYP*', 'SEED')
    # Path to test HMM
    HMM_PATH = os.path.join(TMP_PATH, 'MGYP*', 'HMM')
    # Path to an uncompressed MGnify chunk
    MGNIFY_PATH = os.path.join(ROOT_PATH, 'data', 'mgnify', 'chunk{:03d}.fa.gz')

    # Wrapper for hmmsearch run and evaluation
    def run_hmm_search(hmm_search, hmm_path, chunk_idx=0, cpu=None, z=None):
        # Log search start
        print('Searching for HMM against MGnify...', end='', flush=False)
        # Define path to compressed dataset chunk
        chunk_path = MGNIFY_PATH.format(chunk_idx)
        # Create temporary input file (uncompressed)
        chunk_path = gunzip(
            in_path=chunk_path,  # File to unzip
            out_suffix='.fasta'  # temporary file extension
        )
        # Run hmmsearch using previously created HMM against given MGnify dataset
        results, took = benchmark(
            fn=hmm_search.run,  # Function to run
            hmm_path=hmm_path,  # Path to HMM to search against target dataset
            ds_path=chunk_path,  # Path to target dataset
            domtblout_path='/dev/stdout',  # Path to domtblout
            cpu=cpu,  # How many CPUs to use
            z=z  # Number of sequences in target dataset
        )
        # Debug
        print('DONE in {:.0f} seconds'.format(took))
        # Initialize domtblout results
        matches = list()
        # Loop through every line in table
        for row in iter_domtblout(iter(results.stdout.split('\n'))):
            # Store row
            matches.append(row)
        # Return rows in domtblout
        return matches, took

    # Define a new instance of hmmbuild script
    hmm_build = HMMBuild()
    # Define first seed alignment in list
    seed_path = next(glob.iglob(SEED_PATH))
    # Define a temporary file for HMM
    hmm_path = NamedTemporaryFile(delete=False).name
    # Log build start
    print('Building HMM from SEED alignment...')
    # Run HMM build
    hmm_build.run(
        msa_path=seed_path,
        out_path=hmm_path,
        name='test'
    )

    # Load HMM result from file
    hmm_object = HMM.from_file(hmm_path)
    # Remove temporary hmm file
    os.remove(hmm_path)
    # Debug
    print('Loaded HMM object:')
    print('  name:', hmm_object.name)
    print('  length:', hmm_object.length)
    print('  alphabet:', hmm_object.alphabet)
    print()

    # Benchmarking: multiprocessing through Dask
    # Make new cluster
    cluster = LSFCluster(cores=1, memory='8 GB', walltime='2:00', use_stdin=True)
    # Require cluster resources
    cluster.scale(10)
    # Instantiate new client
    client = Client(cluster)

    # Define a new instance of hmm search script
    hmm_search = HMMSearch()

    # Initialize list of hmm objects
    hmm_lib = list()
    # Debug
    print(HMM_PATH)
    # Loop trhough each available HMM file
    for hmm_path in glob.iglob(HMM_PATH):

        # Load HMM model from file
        hmm_obj = HMM.from_file(hmm_path)
        # Retrieve current HMM model attributes
        hmm_name = hmm_obj.name
        hmm_len = hmm_obj.length

        # Verbose
        print('Loaded HMM model at {:s}...'.format(hmm_path))
        print('name: {:s}'.format(hmm_name))
        print('length: {:d}'.format(hmm_len))

        # Define wethe the
        fit = HMMSearch.get_resources(max_mem=8, model_len=hmm_len)

        # Case it does not fit memory
        if not fit:
            # Verbose
            print('hmmsearch requires more resources than available:', end=' ')
            print('model {:s} has been discarded')
            print()
            # Skip iteration
            continue

        # Verbose
        print('hmmsearch fits available resources:', end=' ')
        print('model {:s} added to library')
        print()

        # Add current model to library
        hmm_lib.append(hmm_path)

    # Define number of HMM models to test together: one, half, all
    num_hmm = [1, len(hmm_lib)//2, len(hmm_lib)]
    # Define number of chunks to be tested together
    num_chunks = [1, 10, 100]
    # Define all possible parameters combinations
    params = product(num_hmm, num_chunks)
    # Looping through each parameter combination
    for max_hmm, max_chunks in params:

        # Initialize timers
        time_beg, time_end = time.time(), None

        # Initialize current HMM library empty file
        hmm_lib_path = NamedTemporaryFile(delete=False).name
        # Fill current HMM library file
        HMM.concat(in_paths=hmm_lib[:max_hmm], out_path=hmm_lib_path)

        # Verbose
        print('Searching {:d} HMM models against {:d} dataset chunks...'.format(
            max_hmm, max_chunks
        ), end='', flush=False)

        # Initialize futures
        futures = list()
        # Loop through given dataset chunks
        for i in range(0, max_chunks):
            # Submit function
            future = client.submit(
                run_hmm_search,  # Run hmmsearch against i-th dataset chunk
                hmm_search=hmm_search,  # Instance of HMMSearch object
                hmm_path=hmm_lib_path,  # Define path to HMM
                chunk_idx=i  # Define i-th dataset chunk
            )
            # Save future
            futures.append(future)

        # Get results
        results = client.gather(futures)
        # Define time took to run
        took = [results[i][1] for i in range(len(results))]
        # Define domtblout
        domtblout = [
            results[i][0][j]
            for i in range(len(results))
            for j in range(len(results[i][0]))
        ]

        # Update timer
        time_end = time.time()
        # Verbose
        print('DONE')
        print('Retrieved {:d} matches'.format(len(domtblout)), end=' ', flush=False)
        print('in {:.0f} seconds'.format(time_end - time_beg), end=' ', flush=False)
        print('({:.0f} seconds/dataset)'.format(np.mean(took)), end=' ', flush=False)
        print()

    # print('Estimated amount of memory required by hmmsearch', end=' ')
    # print('is of {:d} Gb over {:d} processes'.format(req_memory, req_cpus))
    #
    # # Define a new instance of hmm search script
    # hmm_search = HMMSearch()
    # # Test hmmsearch on first chunk of dataset
    # results, took = run_hmm_search(hmm_search, hmm_path=HMM_PATH, chunk_idx=0)
    # # Debug
    # print('Retrieved {:d} matches in {:.0f} seconds'.format(
    #     len(results), took
    # ))
    #
    # # Benchmarking: multiprocessing through Dask
    # # Make new cluster
    # cluster = LSFCluster(cores=1, ncpus=1, memory='16 GB', walltime='2:00', use_stdin=True)
    # # Adapt cluster resources
    # cluster.adapt(minimum=1, maximum=100)
    # # Instantiate new client
    # client = Client(cluster)
    #
    # # Debug
    # print('Cluster:', cluster)
    # print('Client:', client)
    #
    # # Try executing parallel HMMs
    # try:
    #     # Loop through number of cpus
    #     for j in range(0, 1):
    #         # Initialize timers
    #         time_beg, time_end = time.time(), None
    #
    #         # Initialize futures
    #         futures = list()
    #         # Loop through 10 dataset chunks
    #         for i in range(0, 100):
    #             # Submit function
    #             future = client.submit(
    #                 run_hmm_search,  # Run hmmsearch against i-th dataset chunk
    #                 hmm_search=hmm_search,  # Instance of HMMSearch object
    #                 hmm_path=HMM_PATH,  # Define path to HMM
    #                 chunk_idx=i,  # Define i-th dataset chunk
    #                 # cpu=j  # Define number of cores
    #             )
    #             # Save future
    #             futures.append(future)
    #
    #         # Get results
    #         results = client.gather(futures)
    #         # Define time took to run
    #         took = [results[i][1] for i in range(len(results))]
    #         # Define domtblout
    #         domtblout = [
    #             results[i][0][j]
    #             for i in range(len(results))
    #             for j in range(len(results[i][0]))
    #         ]
    #
    #         # Update timer
    #         time_end = time.time()
    #         # Debug
    #         print('Retrieved {:d} matches'.format(len(domtblout)), end=' ', flush=False)
    #         print('in {:.0f} seconds'.format(time_end - time_beg), end=' ', flush=False)
    #         print('(mean time is {:.0f} seconds/dataset)'.format(np.mean(took)), end=' ', flush=False)
    #         print('with {:d} CPUs'.format(j), flush=False)
    #
    # # Intercept subprocess error
    # except CalledProcessError as err:
    #     # Show error
    #     print('Subprocess exited with code', err.returncode)
    #     print('after executing', err.cmd)
    #     # Print traceback
    #     traceback.print_exc()
