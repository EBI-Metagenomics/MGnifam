"""
Handle multiple sequence alignments
"""

# Dependencies
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import tempfile
import sys
import os
import re


class MSA(object):

    # Static attributes
    # List of possible residues (static)
    res = np.array([
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
    ])
    # Gap character
    gap = '-'

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
            msa.beg += self.beg + np.sum(self.aln[:, :i] != self.gap, axis=1)
            msa.end += self.end - np.sum(self.aln[:, j+1:] != self.gap, axis=1)
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

    # Load from aligned file
    def from_aln(self, in_path):
        # Initialize entries dictionary (key is accession)
        entries = dict()
        # Read file
        with open(in_path, 'r') as in_file:
            # Loop through each line in file
            for line in in_file:
                # Get current (last) key
                acc = [*entries.keys()][-1] if entries else ''
                # Check if current line is header
                is_header = re.search(r'^([a-zA-Z0-9]+)/(\d+)-(\d+)[ ]*', line)
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
                # Change dots to minus characters
                line = re.sub(r'[\.\-]', '-', line)
                # Retrieve residues: remove all non-letter characters
                res = list(re.sub(r'[^\-a-zA-Z]+', '', line))
                # Append residues to last entry
                entries[acc]['res'] += res
                # print(entries[acc])
        # Get accession (keys)
        acc = [k for k in entries.keys() if k != '']
        # Get residues matrix
        res = list([entries[k]['res'] for k in acc])
        # Store current data (cast to numpy array)
        self.acc = np.array(acc, dtype=np.unicode_)
        self.aln = np.array(res, dtype=np.unicode_)
        self.beg = np.array([entries[k]['beg'] for k in acc], dtype=np.int)
        self.end = np.array([entries[k]['end'] for k in acc], dtype=np.int)
        # Return self, allow chaining
        return self

    # Print out alignment to fasta file
    def to_aln(self, out_path):
        # Get description
        desc = np.array([
            # Set descriptio line
            str(self.acc[i]) + '/' + str(self.beg[i]) + '-' + str(self.end[i])
            # For every line in the alignment
            for i in range(self.aln.shape[0])
        ])
        # Get max len of description prefix
        desc_len = len(max(desc, key=len))
        # Open file buffer
        with open(out_path, 'w') as out_file:
            # Loop through each line in the alignment
            for i in range(self.aln.shape[0]):
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
        # Return self, allows chaining
        return self

    # Plot alignment matrix as heatmap
    def plot_heatmap(self, ax, score='gap', cmap=None):
        # Get msa shape (n rows, m positions/columns)
        n, m = self.aln.shape
        # Initialize the alignment matrix
        aln = self.aln
        # Get colouring
        if score == 'gap':  # Default: search for gaps
            # Make new alignment matrix
            aln = (aln != self.gap)
            # Set default color map
            cmap = 'binary_r' if cmap is None else cmap
        elif score == 'occupancy':  # Color columns wrt occupancy
            # Compute occupancy
            occupancy, _ = MSA.occupancy(aln)
            # Make new alignment matrix
            aln = np.tile(occupancy, (n, 1))
            # Set default color map
            cmap = 'viridis' if cmap is None else cmap
        elif score == 'conservation':  # Color columns wrt consrevation
            # Compute conservation
            conservation, _ = MSA.conservation(aln)
            # Make new alignment matrix
            aln = np.tile(conservation, (n, 1))
            # Set default color map
            cmap = 'viridis' if cmap is None else cmap
        else:  # Case color is not valid
            raise NotImplementedError('Scoring method is not defined!')
        # Plot image
        im = ax.imshow(aln, cmap=cmap)
        # Add colorbar
        _ = plt.gcf().colorbar(im, ax=ax)
        # Make a mask for gaps
        gaps = (self.aln == self.gap)
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
        # Define y axis according to chosen score
        if score == 'occupancy':
            # Compute score as occupancy
            y, _ = MSA.occupancy(self.aln)
        elif score == 'conservation':
            # Compute score as conservation
            y, _ = MSA.conservation(self.aln)
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
        # Define y axis according to chosen score
        if score == 'occupancy':
            # Compute score as occupancy
            y, _ = MSA.occupancy(self.aln)
        elif score == 'conservation':
            # Compute score as conservation
            y, _ = MSA.conservation(self.aln)
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
            encoded[:, :, z] = (x == res[z])
        # Return encoded matrix as integer matrix
        return encoded.astype(np.int)

    @classmethod
    def consensus(cls, x, axis=0):
        # Get one hot encoded matrix (copy, include gaps)
        x = cls.one_hot_encode(x, include_gap=True)
        # Get shapes (n rows, m columns, k features)
        n, m, k = x.shape
        # Define consensus matrix (m alignment columns x 20 residues)
        freqs = np.sum(x, axis=axis) / x.shape[axis]
        # Return consensus matrix
        return freqs

    @classmethod
    def occupancy(cls, x, axis=0):
        """
        x:      numpy matrix (n rows, m columns) on which occupancy must be
                computed as the sum of occupied cells over all the occupied
                positions. A cell is occupied when symbol is not gap (-);
        """
        # Define consensus matrix (m alignment columns x 20 residues)
        freqs = cls.consensus(x, axis=axis)
        # Define occupancy array (sum all consensus without gaps)
        avg = np.sum(freqs[:, :-1], axis=1)
        # Return either occupancy and consensus
        return avg, freqs

    @classmethod
    def conservation(cls, x, axis=0):
        """
        x:      numpy matrix (n-dimensional) on which conservation is computed,
                as the shannon entropy of each considered cell;
        base:   logarithm base used in Shannon entropy computation;
        """
        # Define consensus matrix (m alignment columns x 20 residues)
        freqs = cls.consensus(x, axis=axis)
        # Compute shannon entropy
        entropy = np.nansum(-(freqs[:, :-1] * np.log(freqs[:, :-1])), axis=1)
        # Return entropy for every column
        return entropy, freqs

    @classmethod
    def prettiness(cls, x):
        """
        Prettiness score is the mean conservation multiplied over the minimum
        between number of rows and number of columns.

        Args:
        x:          numpy matrix (n rows, m columns)

        Return:
        (float):    prettiness score
        """
        # Get alignment shape as n rows and m columns
        n, m = x.shape
        # Compute conservation
        csv, _ = cls.conservation(x)
        # Return prettiness
        return np.mean(csv) * min(n, m)

    @classmethod
    def redundancy(cls, x):
        """
        Computes redundancy among each pair of sequence in multiple sequence
        alignment, i.e. the number of total equal rediues with respect to
        the total number of (non-gap) residues
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
                seq_i = (x[i, :] != cls.gap).astype(np.int)
                # Non-gap residues in j-th sequence
                seq_j = (x[j, :] != cls.gap).astype(np.int)
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
    def __init__(self, cmd='muscle', env=os.environ.copy()):
        # Store attributes
        self.cmd = cmd
        self.env = env

    def run(self, sequences, verbose=False):
        """
        Takes as input list of n sequences, then creates a multiple sequence
        alignment with n rows (one per sequence) and a nomber of coulums
        m < m1 + m2 with m1 and m2 are the two longest sequences.

        Args
        sequences (list)    List of strings containing fasta entries
        verbose (bool)      Wether to show verbose log or not

        Return
        (msa.MSA)           Instance of multiple sequence alignment
        """
        # Define a new temporary file containing input sequences
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        # Open temporary fasta file
        with open(fasta_file.name, 'w') as ff:
            # Get number of input sequences
            n = len(sequences)
            # Loop through each input entry
            for i in range(n):
                # Write input fasta sequences to temporary fasta file
                ff.write(sequences[i] + '\n')

        # Open new output temporary file
        out_file = tempfile.NamedTemporaryFile(delete=False)
        # Run Muscle multiple sequence alignment
        cmd_out = subprocess.run(
            capture_output=True,  # Capture command output
            encoding='utf-8',  # Set output encoding to string
            check=True,  # Check execution (error code)
            env=self.env,  # Set command environment
            # Define command line  arguments
            args=[
                # self.cmd, '-quiet',
                self.cmd,
                '-in', fasta_file.name,  # Path to input temporary file
                '-out', out_file.name  # Path to output temporary file
            ]
        )

        # Read output temporary file
        with open(out_file.name, 'rb') as of:
            # Debug
            print(out_file.read().decode('utf-8'))

        # Delete temporary files
        os.remove(fasta_file.name)
        os.remove(out_file.name)

# Unit testing
if __name__ == '__main__':

    # # Read multiple sequence alignment
    # msa = MSA().from_aln('data/MGYP001224746368/SEED')
    # # Check number of aligned sequences
    # print('There are %d aligned sequences with %d columns' % msa.aln.shape)
    # # Check first 5 rows
    # for i in range(0, min(5, msa.aln.shape[0])):
    #     # Print sequence description
    #     print('Sequence nr %d: %s from %d to %d' % (
    #         i+1, msa.acc[i], msa.beg[i], msa.end[i]
    #     ))
    #     # Print actual sequence
    #     print(''.join(msa.aln[i, :]))
    #     print()
    #
    # # # Test one hot encoder
    # # ohe = MSA.one_hot_encode(msa.aln, include_gap=True)
    # # # Check output
    # # print('Input shape is {}'.format(msa.aln.shape))
    # # print('Output shape is {}'.format(ohe.shape))
    # # # Check first 5 rows
    # # for i in range(0, min(5, msa.aln.shape[0])):
    # #     # Print a subsequence
    # #     print(msa.aln[i, 100:105])
    # #     # Print encoded subsequence
    # #     print(ohe[i, 100:105, :])
    #
    # # Compute occupancy score
    # occupancy, consensus = MSA.occupancy(msa.aln)
    # # Compute conservation score
    # conservation, consensus = MSA.conservation(msa.aln)
    # # Check results
    # print('Input alignment has shape {}'.format(msa.aln.shape))
    # print('Occupancy scores has shape {}'.format(occupancy.shape))
    # print('Conservation scores has shape {}'.format(conservation.shape))
    # # Initialize table of results
    # print('Column\tBest residue\t\tConservation\tOccupancy')
    # # Results (consensus is a matrix with shape m positions x 20 amino acids)
    # for j in range(0, min(consensus.shape[0], 100)):
    #     # Get index of amino acid with best consensus
    #     i = np.argmax(consensus[j, :-1])
    #     # Get best results for this column
    #     print('{:d}\t{:s} (freq={:.03f})\t\t{:.03f}\t\t{:.03f}'.format(
    #         j+1,  # Current position/column in msa
    #         MSA.res[i],  # Amino acid with highest frequency
    #         consensus[j, i],  # Highest frequency
    #         conservation[j],  # Conservation
    #         occupancy[j]  # Occupancy
    #     ))
    #
    # # Make complete trimming (keep none)
    # trimmed = msa.trim(0, 0)
    # # Test trimming
    # print('Input alignment has shape {}'.format(msa.aln.shape))
    # print('Trimmed alignment has shape {}'.format(trimmed.aln.shape))
    # # Loop through each line in alignment
    # for i in range(trimmed.aln.shape[0]):
    #     print('{}:\tsequence {}\t{}-{}'.format(
    #         i+1,  # Row number
    #         trimmed.acc[i],  # Sequence accession
    #         trimmed.beg[i],  # Sequence begin
    #         trimmed.end[i]  # Sequence end
    #     ))
    #
    # # Make redundancy matrix
    # redundancy = MSA.redundancy(msa.aln)
    # # Test redundancy
    # print('Input alignment has shape {}'.format(msa.aln.shape))
    # print('Redundacny matrix has shape {}'.format(redundancy.shape))
    # print(redundancy)
    #
    # # Plot alignment as heatmap
    # fig, ax = plt.subplots(1, 1, figsize=(30, 15))
    # _ = ax.set_title('Multiple Sequence Alignment')
    # _ = msa.plot_heatmap(score='conservation', ax=ax)
    # _ = plt.show()
    #
    # # Plot alignment as scatterplot
    # fig, ax = plt.subplots(1, 1, figsize=(30, 15))
    # _ = ax.set_title('Multiple Sequence Alignment score')
    # _ = msa.plot_scatter(score='conservation', ax=ax)
    # _ = plt.show()

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
    muscle.run(sequences=sequences, verbose=True)
