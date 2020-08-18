"""
Handle multiple sequence alignments
"""

# Dependencies
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import os
import re

# Custom dependencies
from src.utils import open_file
from src.sequences import fasta_iter


class MSA(object):

    # Static attributes
    # List of possible residues (static)
    res = np.array([
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    ])
    # Gap characters
    gap = '-.'

    # Constructor
    def __init__(self):
        # Object attributes
        self.aln = np.ndarray((0, 0))  # Numpy matrix containing aligned rows
        # Numpy array containing
        self.acc = np.ndarray(0)  # List containing sequences description
        self.beg = np.ndarray(0)  # Start position of the aligned sequence
        self.end = np.ndarray(0)  # End position of the aligned sequence

    def trim(self, i, j):
        """Trim
        Trim method keeps only positions/columns in alignment from i
        to j (included). If i >= j, it means that no column must be kept.
        """
        # Make a new MSA object
        msa = self.__class__()
        # Update MSA values with slice of the current ones
        msa.aln = self.aln[:, i:j]
        msa.acc = self.acc
        # Initialize begin and end positions with zeroes
        msa.beg = np.zeros_like(self.beg, dtype=np.int)
        msa.end = np.zeros_like(self.end, dtype=np.int)
        # Check if there is at least 1 position remaining
        if (i < j):
            # Update begin and end positions
            msa.beg = self.beg + np.sum(np.isin(self.aln[:, :i], list(self.gap), invert=True), axis=1)
            msa.end = self.end - np.sum(np.isin(self.aln[:, j:], list(self.gap), invert=True), axis=1)
        # Return trimmed alignment
        return msa

    def slice(self, indexes):
        """
        Slices the multiple sequence alignment keeping only sequences indexed
        in input array

        Args:
        slice (list):   list of indexes to keep, the others will be cut out

        Return:
        (msa.MSA):      multiple sequence alignment object without non
                        indexed rows
        """
        # Make a new MSA object
        msa = self.__class__()
        # Update MSA values with slice of the current ones
        msa.aln = self.aln[indexes, :]
        msa.acc = self.acc[indexes]
        msa.beg = self.beg[indexes]
        msa.end = self.end[indexes]
        # Return slices multiple sequence alignment object
        return msa

    @property
    def is_empty(self):
        """
        Checks whether there isn't any row or column in the current msa

        Return:
        (bool):     Wether there is at least one cell in multiple sequence
                    matrix (False) or not (True)
        """
        # Get shape (n rows, m cols) for current msa matrix
        n, m = self.aln.shape
        # Check wether there is at least one cell
        return (n * m) <= 0

    @property
    def shape(self):
        # Case alignment is empty
        if self.is_empty:
            # Return zeroes
            return 0, 0
        # Otherwise return alignment shape
        return self.aln.shape

    # Load from aligned file
    @classmethod
    def from_aln(cls, path):
        """Load MSA from .aln file

        Args
        in_path (str)   Path to MSA .aln file (even if gzipped)

        Return
        (MSA)           Loaded MSA object
        """
        # Initialize entries dictionary (key is accession)
        entries = dict()
        # Read file
        with open_file(path, 'r', 'rt') as file:
            # Loop through each line in file
            for line in file:
                # Get current (last) key
                acc = [*entries.keys()][-1] if entries else ''
                # Check if current line is header
                is_header = re.search(r'^(\S+)/(\d+)-(\d+)\s+', line)
                # If header, retrieve it
                if is_header:
                    # Retrieve accession
                    acc = str(is_header.group(1))
                    beg = int(is_header.group(2))
                    end = int(is_header.group(3))
                    # Instantiate new entry
                    entries[acc] = {'acc':acc, 'beg': beg, 'end': end, 'res':list()}
                    # Remove header from line
                    line = re.sub(r'^[a-zA-Z0-9]+/[\d]+-[\d]+[ ]*', '', line)
                # # Change dots to minus characters
                # line = re.sub(r'[\.\-]', '-', line)
                # Retrieve residues: remove all non-letter characters
                res = list(re.sub(r'[^\-\.a-zA-Z]+', '', line))
                # Append residues to last entry
                entries[acc]['res'] += res
                # print(entries[acc])
        # Get accession (keys)
        acc = [k for k in entries.keys() if k != '']
        # Get residues matrix
        res = list([entries[k]['res'] for k in acc])
        # Init new MSA object
        msa = cls()
        # Store current data (cast to numpy array)
        msa.acc = np.array(acc, dtype=np.unicode_)
        msa.aln = np.array(res, dtype=np.unicode_)
        msa.beg = np.array([entries[k]['beg'] for k in acc], dtype=np.int)
        msa.end = np.array([entries[k]['end'] for k in acc], dtype=np.int)
        # Return self, allow chaining
        return msa

    # Print out alignment to fasta file
    def to_aln(self, path):
        # Initialize description column
        desc = list()
        # Initialize description length (maximum number of characters)
        desc_len = 0
        # Loop through all rows in alignment
        for i in range(self.shape[0]):
            # Retrieve parameters for making description string
            acc, beg, end = self.acc[i], self.beg[i], self.end[i]
            # Make description string
            desc_str = r'{:s}/{:d}-{:d}'.format(acc, beg, end)
            # Update maximum description length
            desc_len = max(len(desc_str), desc_len)
            # Store description string
            desc.append(desc_str)

        # Open file buffer
        with open(path, 'w') as out_file:
            # Loop through each line in the alignment
            for i in range(self.shape[0]):
                # Format output line
                line = '{desc:s} {spaces:s}{aln:s}\n'.format(
                    # Description prefix
                    desc=desc[i],
                    # Eventually add some more spaces
                    spaces=' ' * (desc_len - len(desc[i])),
                    # Add alignment string
                    aln=''.join(self.aln[i, :])
                )
                # Write sequence
                out_file.write(line)

    # Load from fasta file
    @classmethod
    def from_fasta(cls, path, regex=r'^>(.*)[\n\r]*$'):
        """ Load alignment from .fasta file

        Args
        path (str)          Path to fasta file holding sequences
        regex (str)         Regex used to extract accession number (or a
                            generic id string) from fasta header lines
                            (by default all the line except trailing '>' char)

        Return
        MSA                 Return new Multiple Sequnece Alignment object
        """
        # Initialize output MSA object
        msa = cls()
        # Initialize alignment matrix as list (of numpy arrays)
        aln = list()
        # Initialize accession number, begin and end positions as lists
        acc, beg, end = list(), list(), list()
        # Open file
        with open(path, 'r') as file:
            # Loop through fasta file lines
            for entry in fasta_iter(file):
                # Split entry in header and residues
                header, residues = tuple(entry.split('\n'))
                # Save accession number / id
                acc.append(re.search(regex, header).group(1))
                # Save residues for current aligned sequence
                aln.append(list(residues))

        # Update attributes
        msa.aln = np.array(aln)
        msa.acc = np.array(acc)
        # Count number of non-gap reisudes for each aligned sequence
        res = np.sum(np.isin(msa.aln, list(cls.gap), invert=True), axis=1)
        # Define begin position of each aligned sequence
        msa.beg = np.ones_like(res, dtype=np.int)
        # Define end position of each aligned sequence
        msa.end = res.astype(np.int)

        # Return reference to current object (allows chaining)
        return msa

    # Load from stockholm file (hmmalign output)
    @classmethod
    def from_sto(cls, path):
        """ Load alignment from .sto stockholm file

        Args
        path (str)      Path to .stockholm file holding alignment

        Return
        MSA             Return new Multiple Sequnece Alignment object
        """
        # Initialize new MSA object
        msa = cls()
        # Initialize list of sequences
        aln = list()
        # Initialize list of accession numbers
        acc = list()
        # Initialize offset list
        beg, end = list(), list()

        # Initialize an intermediate dictionary
        sequences_dict = dict()
        # Open file
        with open(path, 'r') as file:
            # Loop through stockholm file lines
            for line in file:
                # Define expected line format to match
                match = re.search(r'^([^#]\S*)\s+(\S+)', line)
                # Case current line does not match expected format
                if not match:
                    # Go to next iteration
                    continue
                # Get raw accession
                accession = match.group(1)
                # Get residues list
                residues = list(match.group(2))
                # Initialize entry
                sequences_dict.setdefault(accession, residues)
                # Update entry
                sequences_dict[accession] += list(residues)

        # Update list of sequences
        for accession, residues in sequences_dict.items():
            # Split accession in accession number and residues boundaries
            found = re.search(r'(\S+)/(\d+)-(\d+)', accession)
            # Retrieve accession number and residues boundaries
            acc.append(found.group(1))
            beg.append(int(found.group(2)))
            end.append(int(found.group(3)))
            # Save residues
            aln.append(residues)

        # Update output MSA
        msa.aln = np.array(aln)
        msa.acc = np.array(acc)
        msa.beg = np.array(beg, dtype=np.int)
        msa.end = np.array(end, dtype=np.int)
        # Return Multiple Sequence Alignment class
        return msa

    # Plot alignment matrix as heatmap
    def plot_heatmap(self, ax, score='gap', cmap=None):
        # Get msa shape (n rows, m positions/columns)
        n, m = self.aln.shape
        # Initialize the alignment matrix
        aln = self.aln
        # Case scoring is according to gaps
        if score == 'gap':  # Default: search for gaps
            # Make new alignment matrix
            aln = np.isin(aln, list(self.gap), invert=True)
            # Set default color map
            cmap = 'binary_r' if cmap is None else cmap
        # Otherwise
        else:
            # Compute consensus
            cns = consensus(aln)
            # Color columns according to occupancy
            if score == 'occupancy':
                # Compute occupancy
                occ = occupancy(cns)
                # Make new alignment matrix
                aln = np.tile(occ, (n, 1))
                # Set default color map
                cmap = 'viridis' if cmap is None else cmap
             # Color columns according to consrevation
            elif score == 'conservation':
                # Compute conservation
                csv = conservation(cns)
                # Make new alignment matrix
                aln = np.tile(csv, (n, 1))
            else:  # Case color is not valid
                raise NotImplementedError('Scoring method is not defined!')
            # Set default color map
            cmap = 'viridis' if cmap is None else cmap
        # Plot image
        im = ax.imshow(aln, cmap=cmap)
        # Add colorbar
        _ = plt.gcf().colorbar(im, ax=ax)
        # Make a mask for gaps
        gaps = np.isin(self.aln, list(self.gap))
        # Add alpha:
        _ = ax.imshow(
            # Set balck only gap positions
            np.ma.masked_where(gaps == 0, gaps),
            # Set colouring scheme
            cmap='viridis',
            # Set alpha level
            alpha=.5
        )
        # Set x and y axis ticks/labels (same)
        x, y = np.arange(0, m, 50), np.arange(0, n, 10)
        # We want to show all ticks
        _ = ax.set_xticks(x)
        _ = ax.set_yticks(y)
        # Rotate the tick labels and set their alignment.
        _ = plt.setp(ax.get_xticklabels(), rotation=90, ha='right', rotation_mode='anchor')
        # Return updated axis
        return ax

    # Plot alignment score as scatterplot
    def plot_scatter(self, ax, score='occupancy'):
        # Get msa shape (n rows, m positions/columns)
        n, m = self.aln.shape
        # Define x axis (alignment positions)
        x = np.arange(0, m)
        # Compute consensus
        cns = consensus(self.aln)
        # Define y axis according to chosen score
        if score == 'occupancy':
            # Compute score as occupancy
            y = occupancy(cns)
        elif score == 'conservation':
            # Compute score as conservation
            y = conservation(cns)
        else:  # Case no valid score
            raise NotImplementedError('Scoring method is not defined!')
        # Make scatterplot
        _ = ax.plot(x, y, '.-')
        # Set x ticks
        _ = ax.set_xticks(np.arange(0, m, 50))
        # Rotate the tick labels and set their alignment.
        _ = plt.setp(ax.get_xticklabels(), rotation=90, ha='right', rotation_mode='anchor')
        # Return updated axis
        return ax

    # Plot alignment score as histogram
    def plot_hist(self, ax, score='occupancy', bins=100, density=True, **kwargs):
        # Get msa shape (n rows, m positions/columns)
        n, m = self.aln.shape
        # Compute consensus
        cns = consensus(self.aln)
        # Define y axis according to chosen score
        if score == 'occupancy':
            # Compute occupancy score
            y = occupancy(cns)
        elif score == 'conservation':
            # Compute conservation score
            y = conservation(cns)
        else:  # Case no valid score
            raise NotImplementedError('Scoring method is not defined!')
        # Make scatterplot
        _ = ax.hist(y, bins=bins, density=density, **kwargs)
        # Return updated axis
        return ax

    @classmethod
    def one_hot_encode(cls, x, include_gap=True):
        # Get residues array, eventually add gap
        res = np.append(cls.res, cls.gap) if include_gap else cls.res
        # Get number of residues
        k = res.shape[0]
        # Get input matrix shape as n rows and m columns
        n, m = x.shape
        # Initialize output matrix (n rows, m columns, k features)
        encoded = np.zeros((n, m, k), dtype=np.int)
        # Loop through each possible residues
        for z in range(k):
            # Set 1 if values match current residue letter
            encoded[:, :, z] = np.isin(x, res[z])
        # Return encoded matrix as integer matrix
        return encoded.astype(np.int)


# Compute consensus
def consensus(x, axis=0):
    # Get one hot encoded matrix (copy, include gaps)
    x = MSA.one_hot_encode(x, include_gap=True)
    # Get shapes (n rows, m columns, k features)
    n, m, k = x.shape
    # Define consensus matrix (m alignment columns x 20 residues)
    freq = np.sum(x, axis=axis) / x.shape[axis]
    # Return consensus matrix
    return freq


# Compute occupancy
def occupancy(x):
    """Compute occupancy score

    Args
    x (np.ndarray)      Consensus matrix
    axis (int)          Axis along which occupancy score must be computed
                        (0 for rows, 1 for columns)

    Return
    (np.array)          Occupancy score array
    """
    # Define occupancy array (sum all consensus without gaps)
    occ = np.sum(x[:, :-1], axis=1)
    # Return either occupancy and consensus
    return occ


# Compute conservation
def conservation(x):
    """Compute conservation bit-score

    Args
    x (np.ndarray)      Consensus matrix
    axis (int)          Axis along which conservation bit-score must be computed
                        (0 for rows, 1 for columns)

    Return
    (np.array)          Conservation bit-score array
    """
    # Compute shannon entropy
    entropy = np.nansum(-(x[:, :-1] * np.log(x[:, :-1])), axis=1)
    # Return entropy for every column
    return entropy


# Compute prettiness
def prettiness(x, n, m):
    """Compute prettiness score
    Prettiness score is the mean conservation multiplied over the minimum
    between number of rows and number of columns.

    Args
    x (np.ndarray)  Conservation array for a specific MSA
    n (int)         Number of rows in MSA
    m (int)         Number of columns in MSA

    Return:
    (float)         prettiness score
    """
    # Return prettiness
    return np.mean(x) * min(n, m)


# Compute redundancy
def redundancy(x):
    """Compute pairwise redundancy

    Compute pairwise redundancy among each pair of sequences in multiple
    sequence alignment, i.e. the number of total equal rediues with respect to
    the total number of (non-gap) residues.

    Args
    x (np.ndarray)      2D matrix representing multiple sequence alignment

    Return
    (np.array)          Squared matrix containing pairwise redundancy scores
    """
    # Get shape of input alignment matrix (n rows, m columns)
    n, m = x.shape
    # Make a new matrix with shape n x n
    redundancy = np.identity(n, dtype=np.float)
    # Loop through each row in redundancy matrix
    for i in range(0, n):
        # Loop through each column in redundancy matrix
        for j in range(i+1, n):
            # Non-gap residues in i-th sequence
            seq_i = (x[i, :] != MSA.gap).astype(np.int)
            # Non-gap residues in j-th sequence
            seq_j = (x[j, :] != MSA.gap).astype(np.int)
            # Non-gap and equal residues among i-th and j-th sequence
            seq_i_j = seq_i * seq_j * (x[i, :] == x[j, :]).astype(np.int)
            # Equal residues wrt number of residues in i-th seqeunce
            redundancy[i, j] = np.sum(seq_i_j) / np.sum(seq_i)
            # Equal residues wrt number of residues in j-th seqeunce
            redundancy[j, i] = np.sum(seq_i_j) / np.sum(seq_j)
    # Return redundancy matrix
    return redundancy


class Muscle(object):

    # Constructor
    def __init__(self, cmd=['muscle'], env=os.environ.copy()):
        # Store attributes
        self.cmd = cmd
        self.env = env

    # Just a wrapper for run function
    def __call__(self, *args, **kwargs):
        # Just call run method
        return self.run(*args, **kwargs)

    # def run(self, sequences, acc_regex='^>(\S+)', verbose=False):
    #     """
    #     Takes as input list of n sequences, then creates a multiple sequence
    #     alignment with n rows (one per sequence) and a nomber of coulums
    #     m < m1 + m2 with m1 and m2 are the two longest sequences.
    #
    #     Args
    #     sequences (list)    List of strings containing fasta entries
    #     acc_regex (str)     Regex used to extract accession number (or a
    #                         generic id string) from fasta header lines
    #                         (by default first string after '>' character)
    #     verbose (bool)      Wether to show verbose log or not
    #
    #     Return
    #     (msa.MSA)           Instance of multiple sequence alignment
    #     """
    #     # Define a new temporary file containing input sequences
    #     fasta_file = tempfile.NamedTemporaryFile(delete=False)
    #     # Open temporary fasta file
    #     with open(fasta_file.name, 'w') as ff:
    #         # Get number of input sequences
    #         n = len(sequences)
    #         # Loop through each input entry
    #         for i in range(n):
    #             # Write input fasta sequences to temporary fasta file
    #             ff.write(sequences[i] + '\n')
    def run(self, fasta_path='', out_path='', quiet=True):
        """ Make MUSCLE multiple sequence alignment

        Args
        fasta_path (str)            Path to input fasta file
        out_path (str)              Path where to store output alignment

        Raise
        (CalledProcessError)        In case it was impossible to run MUSCLE
        (FileNotFoundError)         In case input fasta file cen not be found
        """

        # Check input fasta file path
        if not os.path.isfile(fasta_path):
            # Raise exception
            raise FileNotFoundError(' '.join([
                'could not find input fasta file',
                'at {:s}'.format(fasta_path)
            ]))

        # # Open new output temporary file
        # out_file = tempfile.NamedTemporaryFile(delete=False)
        # # Run Muscle multiple sequence alignment
        # cmd_out = subprocess.run(
        #     capture_output=True,  # Capture command output
        #     encoding='utf-8',  # Set output encoding to string
        #     check=True,  # Check execution (error code)
        #     env=self.env,  # Set command environment
        #     # Define command line  arguments
        #     args=[
        #         *self.cmd, '-quiet',  # Do not intercept command line output
        #         # self.cmd,
        #         '-in', fasta_file.name,  # Path to input temporary file
        #         '-out', out_file.name  # Path to output temporary file
        #     ]
        # )

        # Initialize command line arguments
        args = self.cmd
        # Set quiet
        args += ['-quiet'] if quiet else []
        # Add input file path
        args += ['-in', fasta_path] if fasta_path else []
        # Add output file path
        args += ['-out', out_path] if out_path else []

        # # Initialize output multiple sequnece alignment
        # out_msa = MSA()
        # # Read output temporary file
        # with open(out_file.name, 'r') as of:
        #     # # Debug
        #     # print(out_file.read().decode('utf-8'))
        #     # Read multiple sequence alignment from fasta file
        #     out_msa = out_msa.from_fasta(of, acc_regex=acc_regex)
        #
        # # Delete temporary files
        # os.remove(fasta_file.name)
        # os.remove(out_file.name)
        #
        # # Return output multiple sequence alignment
        # return out_msa

        # Run command line script
        return subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # STDOUT/STDERR encoding
            check=True,  # Check that script worked properly
            env=self.env,  # Environmental variables
            args=args
        )

    @staticmethod
    def parse(self, out_path):
        """ Load file as MSA instance

        Args
        out_path (str)          Path to Muscle output file

        Raise
        (FileNotFoundError)     In case output file is not valid
        """
        raise NotImplementedError


# Unit testing
if __name__ == '__main__':

    # Define root path
    ROOT_PATH = os.path.dirname(__file__) + '/../..'
    # Define examples path
    EXAMPLES_PATH = ROOT_PATH + '/tmp/examples/MGYP*'

    # Load all examples paths
    examples_path = glob(EXAMPLES_PATH)

    # Define alignment path
    aln_path = os.path.join(examples_path[0], 'SEED')
    # Read multiple sequence alignment
    msa = MSA.from_aln(aln_path)
    # Check number of aligned sequences
    print('There are {:d} aligned sequences'.format(msa.aln.shape[0]), end=' ')
    print('with {:d} columns'.format(msa.aln.shape[1]), end=' ')
    print('in alignment {:s}'.format(aln_path))
    # Check first 5 rows
    for i in range(0, min(5, msa.aln.shape[0])):
        # Print sequence description
        print('Sequence nr %d: %s from %d to %d' % (
            i+1, msa.acc[i], msa.beg[i], msa.end[i]
        ))
        # Print actual sequence
        print(''.join(msa.aln[i, :]))
        print()

    # # Define alignment path
    # sto_path = os.path.join(examples_path[0], 'ALIGN.sto')
    # # Open alignment file
    # with open(sto_path, 'r') as sto_file:
    #     # Read multiple sequence alignment from stockholm
    #     msa = MSA.from_sto(sto_file)
    #     # Check number of aligned sequences
    #     print('There are {:d} aligned sequences'.format(msa.aln.shape[0]), end=' ')
    #     print('with {:d} columns'.format(msa.aln.shape[1]), end=' ')
    #     print('in alignment {:s}'.format(sto_path))
    #     # Check first 5 rows
    #     for i in range(0, min(5, msa.aln.shape[0])):
    #         # Print sequence description
    #         print('Sequence nr %d: %s from %d to %d' % (
    #             i+1, msa.acc[i], msa.beg[i], msa.end[i]
    #         ))
    #         # Print actual sequence
    #         print(''.join(msa.aln[i, :]))
    #         print()

    # # Test one hot encoder
    # ohe = MSA.one_hot_encode(msa.aln, include_gap=True)
    # # Check output
    # print('Input shape is {}'.format(msa.aln.shape))
    # print('Output shape is {}'.format(ohe.shape))
    # # Check first 5 rows
    # for i in range(0, min(5, msa.aln.shape[0])):
    #     # Print a subsequence
    #     print(msa.aln[i, 100:105])
    #     # Print encoded subsequence
    #     print(ohe[i, 100:105, :])

    # Compute consensus
    cns = consensus(msa.aln)
    # Compute occupancy score
    occ = occupancy(cns)
    # Compute conservation score
    csv = conservation(cns)
    # Check results
    print('Input alignment has shape {}'.format(msa.aln.shape))
    print('Occupancy scores has shape {}'.format(occ.shape))
    print('Conservation scores has shape {}'.format(csv.shape))
    # Initialize table of results
    print('Column\tBest residue\t\tConservation\tOccupancy')
    # Results (consensus is a matrix with shape m positions x 20 amino acids)
    for j in range(0, min(cns.shape[0], 100)):
        # Get index of amino acid with best consensus
        i = np.argmax(cns[j, :-1])
        # Get best results for this column
        print('{:d}\t{:s} (freq={:.03f})\t\t{:.03f}\t\t{:.03f}'.format(
            j+1,  # Current position/column in msa
            MSA.res[i],  # Amino acid with highest frequency
            cns[j, i],  # Highest frequency
            csv[j],  # Conservation
            occ[j]  # Occupancy
        ))

    # Make complete trimming (keep none)
    trimmed = msa.trim(10, 50)
    # Test trimming
    print('Input alignment has shape {}'.format(msa.aln.shape))
    print('Trimmed alignment has shape {}'.format(trimmed.aln.shape))
    # Loop through each line in alignment
    for i in range(trimmed.aln.shape[0]):
        print('{}:\tsequence {}\t{}-{}'.format(
            i+1,  # Row number
            trimmed.acc[i],  # Sequence accession
            trimmed.beg[i],  # Sequence begin
            trimmed.end[i]  # Sequence end
        ))

    # # Make redundancy matrix
    # redundancy = MSA.redundancy(msa.aln)
    # # Test redundancy
    # print('Input alignment has shape {}'.format(msa.aln.shape))
    # print('Redundacny matrix has shape {}'.format(redundancy.shape))
    # print(redundancy)

    #

    # Plot alignment as heatmap
    fig, ax = plt.subplots(1, 1, figsize=(30, 15))
    _ = ax.set_title('Multiple Sequence Alignment')
    _ = msa.plot_heatmap(score='occupancy', ax=ax)
    _ = plt.show()

    # Plot alignment as scatterplot
    fig, ax = plt.subplots(1, 1, figsize=(30, 15))
    _ = ax.set_title('Multiple Sequence Alignment score')
    _ = msa.plot_scatter(score='occupancy', ax=ax)
    _ = plt.show()

    # Define multiple sequence alignment input test sequences
    sequences = [
        ">MGYP001026457211 PL=10 UP=0 BIOMES=0000000011000 LEN=188 CR=0\nNTTCENADRLGTIPADSYAKVNSPSFLGIPLTTTPKPDAPSTQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWESVFSEVALSINALNAAVDQPTSPNDYSGQLKFTGKRKITALNLTNIKTTTSEYATVIGMRADNKELAYEFICIDNYIYMRTGKGDTWNNAISIISKFTSF",
        ">MGYP000848664103 PL=00 UP=0 BIOMES=0000000011000 LEN=297 CR=1\nMAEYSSELDKITYAELALSLQNTIKNNLAHTKDQVIHVTQEDKNKWNQISDIPEATETKKGALTPQEKIKLKNIEERANNYTHPTSGVTAGQYIQVEVNAEGHVVAGHNPTKINTTCENADRLGTIPADSYAKVNSPSFLGIPLTTTPKPDAPSTQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWQSVFSEVALSINALNAAVDQPTSPNDYSGQLKFTGKRKITALKLTNIKATTSEYATVIGMRADNRELAYEFICIDNYIYMRTGKGDTWNNAISIIKD",
        ">MGYP001028524128 PL=01 UP=0 BIOMES=0000000011000 LEN=26 CR=0\nMAEYSSELDKITYAELALSLQNTIKN",
        ">MGYP000854299633 PL=00 UP=0 BIOMES=0000000011000 LEN=297 CR=0\nMAEYSSELDKITYAELALSLQNTIKSNLAHTKDQVIHVTQEDKNKWNQISDIPEATETKKGALTPQEKIKLKNIEERANNYTHPTSGVTAGQYIQVEVNAEGHVVAGHNPTKINTTCENADRLGTIPADSYAKVNSPSFLGIPLTPTPKLDAPASQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWQSVFSEVALSINALNTAVDQPTSPNDYSGQFKFTGKRKITALKLTNIKATTSEYATVIGMRADNRELAYEFICIDNYIYMRTGKGDTWNNAISIIKD"
    ]
    # Define a Muscle instance
    muscle = Muscle()
    # Run Muscle with input test sequences
    msa = muscle.run(sequences, acc_regex=r'^>(\S+)')
    # Debug
    print(msa.aln)
    # Check number of aligned sequences
    print('There are %d aligned sequences with %d columns' % tuple(msa.aln.shape))
    # Check all test rows
    for i in range(msa.aln.shape[0]):
        # Print sequence description
        print('Sequence nr %d: %s from %d to %d' % (
            i+1, msa.acc[i], msa.beg[i], msa.end[i]
        ))
        # Print actual sequence
        print(''.join(msa.aln[i, :]))
        print()
