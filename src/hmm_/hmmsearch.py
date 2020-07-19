# Dependencies
import subprocess
import sys
import os
import re

from src.hmm.hmmer import HMMER


class HMMSearch(HMMER):

    # Constructor
    def __init__(self, cmd=['hmmsearch'], env=os.environ.copy()):
        # Call parent constructor
        super(HMMSearch, self).__init__(cmd=cmd, env=env)

    def run(
        self, hmm_path, target_path, out_path='/dev/null', num_cpus=None,
        seq_e=1000, dom_e=None, seq_z=None, dom_z=None,
        tblout_path='', domtblout_path='', pfamtblout_path=''
    ):
        """Run hmmsearch with some options

        Args
        hmm_path (str)          Path to Hidden Markov Model (HMM) previously
                                created through hmmbuild
        target_path (str)       Path to dataset against which given HMM must be
                                searched
        out_path (str)          Path to output file
        num_cpus (int)          Number of cpus to be used in multithreading
        seq_e (float)           Maximum e-value threshold for sequence match
        dom_e (float)           Maximum e-value threshold for domain match
        seq_z (int)             Z value (length of the dataset) to be used
                                in sequences e-value calculation
        dom_z (int)             Z value (length of the dataset) to be used
                                in domains e-value calculation
        cut_ga                  Use the GA (GAthering) threshold, must be set in
                                alignment .sto (Stockholm) output file
        cut_nc
        cut_tc
        tblout_path (str)       Path to tabular per-sequence output
        domtblout_path (str)    Path to tabular per-domain output
        pfamtblout_path (str)   Path to pfam-formatted tabular output

        Return
        (CompletedProcess)      Instance returned by subprocess.run(...) for a
                                successfully completed process (contains stdout)

        Raise
        (CalledProcessError)    Raised in case hmmbuild process fails
        (Error)                 Unknown error
        """
        # Set hmmsearch executable
        cmd = [*self.cmd]

        # Set output path, if any
        cmd += ['-o', out_path] if out_path else []

        # Set number of cpus used, if any
        cmd += ['--cpu', str(num_cpus)] if num_cpus else []

        # Set sequences e-value
        cmd += ['-E', str(seq_e)] if seq_e is not None else []
        # Set domains e-value
        cmd += ['-domE', str(dom_e)] if dom_e is not None else []
        # Set sequences Z value
        cmd += ['-Z', str(seq_z)] if seq_z is not None else []
        # Set domains Z value
        cmd += ['-domZ', str(dom_z)] if dom_z is not None else []

        # Set tabular formats
        cmd += ['--tblout', tblout_path] if tblout_path else []
        cmd += ['--domtblout', domtblout_path] if domtblout_path else []
        cmd += ['--pfamtblout', pfamtblout_path] if pfamtblout_path else []

        # Set input hmm file
        cmd += [hmm_path]
        # Set input dataset path
        cmd += [target_path]

        # Run HMM search
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )

    @staticmethod
    def get_cpus(max_memory, model_len, max_cpus=1, longest_seq=4e04, num_bytes=48):
        """Computes the resources needed to run hmmscan

        Number of cpus is the maximum number of cpus satisfying the following:
        max     num_cpus
        s.t.    (model_len * long_seq * num_bytes * num_cpus) / 1e09 > mem_threshold
        where   mem_threshold is the maximum amount of memory that can be allocated,
                model_len is the length of the model (found in HMM file),
                long_seq is the length of the longest sequence compared with HMM,
                num_bytes is the number of bytes in the dp

        Args
        max_memory (int)        Maximum number of Bytes that can be allocated to
                                run hmmscan (on a single LSF job)
        longest_seq (int)       Length of the longest sequence in hmmscan target
                                dataset
        model_len (int)         Length of the HMM model
        max_cpus (int)          Maximum number of available CPUs
        num_bytes (int)         Number of bytes in DP

        Return
        (int)                   Number of cpus needed (0 if allocated resources
                                are not fitting computational ones)
        """
        # Compute number of Gb required by a single core
        req_memory = math.ceil(model_len * longest_seq * num_bytes)
        # Loop from maximum number of CPUS to minimum (1)
        for num_cpus in range(max_cpus, 0, -1):
            # Check memory allocated to a single core
            core_memory = max_memory // num_cpus
            # Case allocated memory does not fit required
            if core_memory < req_memory:
                continue  # Skip iteration
            # Otherwise, return current number of cpus
            return num_cpus
        # If flow reaches this point, no cpu can be allocated
        return 0


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
