"""HMM scan and HMM build

This library wraps HMMER functions for either building an Hidden Markov Model
(HMM) starting from a Multiple Sequence Alignment (MSA) and for searching a
HMM against a given dataset (e.g. UniProt).
"""


# Dependecies
import subprocess
import tempfile
import json
import sys
import os
import re


# Custom dependencies
import src.dataset as ds
from src.msa import MSA


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
        self, hmm_path, ds_path, out_path='\\dev\\stdout', cpu=None, z=None,
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

        # Open multiple sequence alignment
        

        # Run HMM search
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )

    @staticmethod
    def iter_domtblout(file):

        # Skip first line
        _ = next(file)
        # Save second line (header, contains column names)
        header_line = re.sub(r'[\n\r]', '', next(file))
        header_line = re.sub(r'^#', ' ', header_line)  # Handle first character
        # Save third line (bounds, define wether a column starts and ends)
        bounds_line = re.sub(r'[\n\r]', '', next(file))
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
        for line in file:
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


# TODO domtblout iterator
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
def iter_tblout(iter):
    raise NotImplementedError

# TODO pfamtblout iterator
def iter_pfamtblout(iter):
    raise NotImplementedError


# Unit testing
if __name__ == '__main__':

    # Define project root path
    ROOT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
    # Define path to temporary folder
    TMP_PATH = os.path.join(ROOT_PATH, 'tmp')
    # Path to test seed file
    SEED_PATH = os.path.join(TMP_PATH, 'MGYP001224746368', 'SEED')
    # Path to test HMM
    HMM_PATH = os.path.join(TMP_PATH, 'test.hmm')
    # Path to an uncompressed MGnify chunk
    MGNIFY_PATH = os.path.join(ROOT_PATH, 'data', 'mgnify', 'chunk000.fa.gz')
    # Path to domtblout file
    DOMTBLOUT_PATH = os.path.join(TMP_PATH, 'test.domtblout')

    # Grab an example MSA
    msa = MSA().from_aln(SEED_PATH)
    # Debug
    print('Input MSA has shape: {}'.format(msa.aln.shape))

    # Define a new instance of hmm build script
    hmm_build = HMMBuild()
    # Log build start
    print('Building HMM from SEED alignment...')
    # Run HMM build
    hmm_build.run(
        msa_path=SEED_PATH,
        out_path=HMM_PATH,
        name='test'
    )

    # Define a new instance of hmm search script
    hmm_search = HMMSearch()
    # Log search start
    print('Searching for HMM against MGnify...')
    # Create temporary input file (uncompressed)
    ds_path = ds.gunzip(MGNIFY_PATH)
    # Run hmmsearch using previously created HMM against given MGnify dataset
    hmm_search.run(
        hmm_path=HMM_PATH,
        ds_path=ds_path,
        domtblout_path=DOMTBLOUT_PATH
    )
    # Remove previously created uncompressed temp file
    os.remove(ds_path)

    # DEBUG Iterate through every line in resulting file, get row dict
    print('HMM search results:')
    # Open buffer to hmmsearch output file
    with open(DOMTBLOUT_PATH, 'r') as file:
        # Define row iterator
        row_iter = enumerate(iter_domtblout(file))
        # Go through each row dict
        for i, row, in row_iter:
            # Print current row
            print('{:d}-th row: {}'.format(i+1, row))
