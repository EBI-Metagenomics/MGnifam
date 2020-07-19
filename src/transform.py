# Common dependencies
import matplotlib.pyplot as plt
import numpy as np

# Custom dependecies
from src.msa import MSA
from src.msa import consensus
from src.msa import occupancy
from src.msa import conservation


class Transform(object):
    """
    Takes as input a MSA object and allows to transform it. This operation involves
    the following steps:
    1. Identification of N and C termial, acccording to occpuancy or
    other scores
    2. Cut according to N and C terminal
    3. Cut out sequences whose similarity is above a certain thershold
    """

    # Wrapper for transform
    def __call__(self, *args, **kwargs):
        return self.transform(*args, **kwargs)

    # Apply transformation
    def transform(self, msa):
        """
        Applies transformation to multiple sequence alignment

        Args
        msa (msa.MSA):  Multiple sequence alignment

        Return
        msa (msa.MSA):  Transformed multiple sequence alignment
        """
        raise NotImplementedError

    # Retrieve threshold function
    def get_threshold(self, threshold):
        """
        Takes an input a threshold, returns a function which has as parameter
        the scoring array and returns a single float value (threshold)

        Args:
        threshold (str/function):   Either a pre defined function name or
                                    a custom function, taking as input x
                                    scoring array and returning float threshold

        Return:
        threshold_fn (function):    Function computing the threshold, whose
                                    only parameter is the scoring array
        """
        # Initialize threshold function
        threshold_fn = None
        if isinstance(threshold, (int, float, complex)):  # Number
            threshold_fn = lambda x: threshold
        elif threshold == 'mean':  # Default mean threshold
            threshold_fn = lambda x: np.mean(x)
        elif threshold == 'median':  # Default median threshold
            threshold_fn = lambda x: np.median(x)
        elif callable(threshold):  # Function
            threshold_fn = threshold
        else:  # Case no valid function
            raise NotImplementedError('Given threshold function is not valid!')
        # Return threshold function
        return threshold_fn

    # Make plots for assessing transformation results
    @classmethod
    def plot_results(cls, msa, axs):
        """Plot transformation results
        This function takes as input a MSA (either transformed or not) and
        returns some plots (axes) which allow to assess transformation results
        goodness.

        Args
        msa (msa.MSA)   Multiple sequence alignment instance
        axs (plt.axes)  Axes (flat) list where to make plots, ordered:
                        axs[0] is the msa coloured by occupancy,
                        axs[1] is the occupancy distribution,
                        axs[2] is the occupation per position scatter,
                        axs[3] is the msa coloured by conservation,
                        axs[4] is the conservation distribution,
                        axs[5] is the conservation per position scatter

        Return
        (plt.axes)      Updated axes
        """

        # Occupancy msa
        axs[0].set_title('MSA, coloured by occupancy')
        msa.plot_heatmap(score='occupancy', ax=axs[0])
        # Occupancy distribution
        axs[1].set_title('Occupancy distribution')
        msa.plot_hist(score='occupancy', density=True, ax=axs[1])
        # Occupancy scatter
        axs[2].set_title('Occupancy per aligned position')
        msa.plot_scatter(score='occupancy', ax=axs[2])

        # Conservation msa
        axs[3].set_title('MSA, coloured by conservation')
        msa.plot_heatmap(score='conservation', ax=axs[3])
        # Conservation distribution
        axs[4].set_title('Conservation distribution')
        msa.plot_hist(score='conservation', density=True, ax=axs[4])
        # Conservation scatter
        axs[5].set_title('Conservation per aligned position')
        msa.plot_scatter(score='conservation', ax=axs[5])

        # Return updated axes
        return axs


class Compose(Transform):
    """Compose transformer
    Compose Trimmer is a special trimmer which does not implement any specific
    trim method but uses the ones in the given list of Trimmer objects.
    Trimmer objects in list are applied to given msa sequentially.
    """

    def __init__(self, transforms):
        """Constructor

        Args:
        transforms (list):  List containing trimmer Objects which will be
                            called on given MSA object sequentially once run
        """
        # Save trimmer list
        self.transforms = transforms

    def transform(self, msa):
        """Trim
        Applies transform() method of Transform objects in input list.

        Args:
        msa (msa.MSA):  Input multiple sequence alignment

        Return:
        msa (msa.MSA):  Output multiple sequence alignment, trimmed
        """
        # Get list of Transform objects
        transforms = self.transforms
        # Loop through each Transform object in list
        for i in range(len(transforms)):
            # Define current Trimmer object
            transform = transforms[i]
            # Apply trim method of current Trimmer object
            msa = transform(msa)
        # Return latest multiple sequence alignment
        return msa


class OccupancyTrim(Transform):
    """Occupancy trimming
    Trims multiple sequence alignment according to given occupancy threshold.
    I.e. it sets i (N terminal) and j (C terminal) as the first indexes where
    threshold is exceeded from left and right, respectively.
    """

    def __init__(self, threshold=0.8, inclusive=True):
        """Constructor

        Args:
        threshold (float):  Occupancy threshold for determining N- and C- term
        inclusive (bool):   Wether to include or not the threshold value when
                            checking for N and C terminal
        """
        self.threshold = self.get_threshold(threshold)
        self.inclusive = inclusive

    def transform(self, msa):
        # Get attributes
        threshold = self.threshold
        inclusive = self.inclusive
        # Compute consensus
        cns = consensus(msa.aln)
        # Compute occupancy
        occ = occupancy(cns)
        # Check which values exceed the threshold
        pos = (occ >= threshold(occ)) if inclusive else (occ > threshold(occ))
        # Get indexes where values exceed the threshold
        pos = np.argwhere(pos == 1).flatten()
        # Get N- and C-terminal as first and last indexes, if any
        i, j = (pos[0], pos[-1]+1) if len(pos) > 1 else (0, 0)
        # Return trimmed MSA
        return msa.trim(i, j)


class ConservationTrim(Transform):
    """Conservation trimming
    Trims multiple sequence alignment according to given conservation threshold.
    I.e. it sets i (N terminal) and j (C terminal) as the first indexes where
    threshold is exceeded from left and right, respectively.
    """

    def __init__(self, threshold=0.8, inclusive=True):
        """Constructor

        Args:
        threshold (float):  Conservation threshold for determining N and C terminal
        inclusive (bool):   Wether to include or not the threshold value when
                            checking for N and C terminal
        """
        self.threshold = self.get_threshold(threshold)
        self.inclusive = inclusive

    def transform(self, msa):
        # Get attributes
        threshold = self.threshold
        inclusive = self.inclusive
        # Compute consensus
        cns = consensus(msa.aln)
        # Compute conservation
        csv = conservation(cns)
        # Check which values exceed the threshold
        pos = (csv >= threshold(csv)) if inclusive else (csv > threshold(csv))
        # Get indexes where values exceed the threshold
        pos = np.argwhere(pos == 1).flatten()
        # Get N- and C-terminal as first and last indexes, if any
        i, j = pos[0], pos[-1]+1 if len(pos) else (0, 0)
        # Return trimmed MSA
        return msa.trim(i, j)


class OccupancyFilter(Transform):
    """Occupancy filtering
    Removes from multiple sequence alignment according those sequence which do
    not fit minimal occupancy threshold.
    """

    def __init__(self, threshold=0.5, inclusive=True):
        """Constructor

        Args:
        threshold (float):  Occupancy threshold for input sequences
        inclusive (bool):   Wether to include or not the threshold value when
                            filtering out rows (keep major or equal vs keep
                            only major)
        """
        self.threshold = self.get_threshold(threshold)
        self.inclusive = inclusive

    def transform(self, msa):
        # Get attributes
        threshold = self.threshold
        inclusive = self.inclusive
        # Compute consensus
        cns = consensus(msa.aln, axis=1)
        # Compute occupancy
        occ = occupancy(cns)
        # Check which values exceed the threshold
        pos = (occ >= threshold(occ)) if inclusive else (occ > threshold(occ))
        # Get indexes where values exceed the threshold
        pos = np.argwhere(pos == 1).flatten().tolist()
        # Return subset of given MSA
        return msa.slice(pos)


# class MakeNonRedundant(Transform):
#     """Make non redundant
#     Make sequence alignment non redundant up to a given threshold: checks
#     wether a bigger sequence contains a smaller, redundant sequence whose
#     redundancy score exceeds a given threshold and removes the samller one
#     in that case
#     """
#
#     def __init__(self, threshold=0.8, inclusive=True):
#         """Constructor
#
#         Args
#         threshold (float):  Redundancy threshold for determining wether a
#                             sequence must be kept or not
#         inclusive (bool):   Wether to include or not the threshold value when
#                             checking for redundancy
#         """
#         self.threshold = self.get_threshold(threshold)
#         self.inclusive = inclusive
#
#     def transform(self, msa):
#         """Make non redundant
#         Starts from longest sequences with highest mean redundancy, drop
#         sequences which have redundancy score exceeding the given threshold
#
#         Args
#         msa (msa.MSA):  Multiple sequence alignment
#
#         Return
#         msa (msa.MSA):  Transformed multiple sequence alignment
#         """
#         # Retrieve transformer attributes
#         threshold = self.threshold
#         inclusive = self.inclusive
#         # Compute redundancy matrix
#         red_mat = MSA.redundancy(msa.aln)
#         if inclusive:  # Apply redundancy threshold (inclusive)
#             red_mat = (red_mat >= threshold(red_mat))
#         else:  # Apply redundancy threshold (exclusive)
#             red_mat = (red_mat > threshold(red_mat))
#         # Get number of residues per sequence
#         num_res = np.sum((msa.aln != msa.gap), axis=1)
#         # Get set of sequence all indices available, sorted by length
#         seq_all = set(np.argsort(num_res).tolist())
#         # Initialize empty set of deleted (redundant) sequences
#         seq_del = set()
#         # Loop through all aligned sequences (column) from longest to smallest
#         for j in seq_all:
#             # Case current sequence has been previously discarded
#             if j in seq_del:
#                 continue  # Skip iteration
#             # Get all sequences contained in current sequence
#             to_del = set(np.argwhere(red_mat[:, j]).flatten().tolist())
#             # Remove current sequence index
#             to_del = to_del - set([j])
#             # Add contained sequence indexes to set of discarded ones
#             seq_del = seq_del | to_del
#         # Slice input MSA according to remaining sequences indexes
#         return msa.slice(list(seq_all - seq_del))


# Test
if __name__ == '__main__':

    # Define path to input seed
    SEED_PATH = './tmp/examples/MGYP000050665084/SEED'

    # Define transformation pipeline
    transform = Compose([
        # Exclude regions outside N- and C- terminal
        OccupancyTrim(threshold=0.4, inclusive=True),
        # Exclude sequences with less than half occupancy
        OccupancyFilter(threshold=0.5, inclusive=True)
        # # Make non redundant with 80 percent threshold
        # MakeNonRedundant(threshold=0.8, inclusive=False)
    ])

    # Read test multiple sequence alignment
    pre_trim = MSA.from_aln(SEED_PATH)
    # Debug
    print('Input multiple sequence alignment has shape', pre_trim.aln.shape)

    # Apply non redundancy to input msa
    post_trim = transform(pre_trim)
    # Debug
    print('Filtered multiple sequence alignment has shape', post_trim.aln.shape)

    # Initialize new plot for comparing pre- transformation scores
    fig = plt.figure(constrained_layout=True, figsize=(20, 10))
    # Define axes grid
    grid = fig.add_gridspec(2, 4)
    # Add main title
    fig.suptitle('Pre-trimming multiple sequence alignment (MSA)')
    # Initialize axis
    axs = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2:]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2:])
    ]
    # Make plots in axes
    Transform.plot_results(msa=pre_trim, axs=axs)
    # Show plot
    plt.show()
    plt.close()

    # Initialize new plot for comparing pre- transformation scores
    fig = plt.figure(constrained_layout=True, figsize=(20, 10))
    # Define axes grid
    grid = fig.add_gridspec(2, 4)
    # Add main title
    fig.suptitle('Post-trimming multiple sequence alignment (MSA)')
    # Initialize axis
    axs = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2:]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2:])
    ]
    # Make plots in axes
    Transform.plot_results(msa=post_trim, axs=axs)
    # Show plot
    plt.show()
    plt.close()

    # Debug
    print(SEED_PATH + '_trimmed')

    # Save to disk
    post_trim.to_aln(SEED_PATH + '_trimmed')
    # Reload from disk
    post_trim = MSA.from_aln(SEED_PATH + '_trimmed')

    # Initialize new plot for comparing pre- transformation scores
    fig = plt.figure(constrained_layout=True, figsize=(20, 10))
    # Define axes grid
    grid = fig.add_gridspec(2, 4)
    # Add main title
    fig.suptitle('Post-trimming multiple sequence alignment (MSA), loaded from disk')
    # Initialize axis
    axs = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2:]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2:])
    ]
    # Make plots in axes
    Transform.plot_results(msa=post_trim, axs=axs)
    # Show plot
    plt.show()
    plt.close()
