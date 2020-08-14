# Dependencies
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
from glob import iglob
from time import time
import os
import re


class HMMER(object):

    def __init__(self, cmd, env=os.environ.copy()):
        self.cmd = cmd
        self.env = env

    def run(self, *args, **kwargs):
        raise NotImplementedError

    # Just a wrapper for run(...) method
    def __call__(self, *args, **kwargs):
        return self.run(*args, **kwargs)


class Tblout(object):

    # Define total number of columns
    width = 19

    # Define columns names
    columns = {
        'target_name': 0,
        'query_name': 2,
        'e_value': 4,
        'bit_score': 5
    }

    # Constructor
    def __init__(self, path):
        # Save path to file
        self.path = path
        # Load file from path
        with open(self.path, 'r') as file:
            # Define file iterator
            iter = self.iterator(file)
            # Initialize inner table
            self.table = list()
            # Iterate through file and store table
            for row in iter:
                # Merge by whitespace al columns after current width
                row = row[:self.width - 1] + [' '.join(row[self.width:])]
                # Add new line to inner table
                self.table.append({
                    col: row[i] for col, i in self.columns.items()
                })

    @staticmethod
    def iterator(iterable):
        """ Retrieve table rows as lists of tuples

        Given an input iterable (e.g. a file or a list of strings) representing
        a tabular output file content, retrieves rows as a list: each item in
        the retrieved list is a row, represented by another list whose elements
        are row attributes, formatted in tuples as (non unique key, value).

        Args
        iterable (iterable)         An iterable object holding hmmsearch tabular
                                    output, such as a file instance or a list of
                                    strings

        Return
        (generator)                 List of rows: each row is a list whose
                                    attributes are cell values
        """
        # Loop through each iterable lines
        for i, line in enumerate(iterable):
            # Check if current line is comment
            is_comment = re.search(r'^#', line)
            # Case current line is comment
            if is_comment:
                # Skip line, go to next one
                continue
            # Clean line from strange characters
            line = re.sub(r'[\n\r]$', '', line)
            # Yield line splitted according to whitespaces
            yield re.split(r'\s+', line)

    # Return dictionary associating each query name with a tuple (TC, NC, GA)
    def get_scores(self, e_value=0.01, min_bits=25.0, max_bits=999999.99):
        """ Compute (TC, NC, GA) threhsold bit scores

        If row e-value is below given e-value threshold, then check if
        previous TC score associated with that query name is above current bit
        score: if so, set the latter as new TC for that query name.

        Otherwise, if row e-value is greater or equal than given e-value
        threshold, check if previous NC score for that query name is below
        current bit score: if so, set the latter as new NC for that query name

        Gathering threshold (GA) is computed at each iteration as mean point
        between NC and TC.

        Notice that both TC and NC are lower bounded bi a bit score of 25.0.
        Moreover, NC is ensured to be lower or equal than TC, it cannot be
        greater. For every query name in given results table file, a tuple
        (TC, NC, GA) is made and returned, there won't be query names with no
        associated threshold bit scores, if found in at least one row.

        Args
        e_value (float)     Threshold e-value used to define a significant hit
                            (hits below that value are considered significant)
        min_bits (float)    Minimum bit-score allowed
        max_bits (float)    Maximum bit-score allowed

        Return
        (dict)              Dictionary mapping cluster name to tuple(TC, NC, GA)

        Raise
        (ValueError)        In case there was an error parsing table rows
        """

        # Initialize scores dictionary mapping query name to tuple(TC, NC, GA)
        scores = dict()
        # Loop through each row in current file
        for i, row in enumerate(self.table):

            # Retrieve query name
            curr_qname = str(row['query_name'])
            # Retrieve e_value from current row
            curr_eval = float(row['e_value'])
            # Retrieve bit score from current row
            curr_bits = float(row['bit_score'])

            # Check if current e-value satisfies threshold
            if curr_eval < e_value:

                # Ensure that entry for current query name is set
                # TC is set extremely high (should include every score)
                # NC is set to lowest possible bit score
                # GA will be recomputed at the first iteration
                tc, nc, ga = scores.setdefault(curr_qname, (min_bits, max_bits, min_bits))

                # Case current bit score is greater than current TC
                if curr_bits > tc:
                    # Update current TC
                    tc = curr_bits

                # Case current bit score is lower than current NC
                elif curr_bits < nc:
                    # Update current NC
                    nc = curr_bits

                # Ensure that NC is lower than TC
                nc = min(nc, tc)
                # Compute gathering threshold (GA): midpoint between NC and TC
                # Ensure that it is not greater than 25.0 (max value)
                ga = min(25.0, nc + ((tc - nc) / 2))

                # Update scores for current query name
                scores[curr_qname] = (tc, nc, ga)

        # Return scores dict
        return scores

    # Return dictionary associating query names with significant target names
    def get_hits(self, scores):
        """ Retrieve sequences that meets scores requirements

        Loop through each row in inner tabular output file: if its bit score
        is greater or equal than the GA set for the associated query name,
        then retrieve target name as an element associated to given query name.

        Args
        scores (dict)       Dictionary associating to each query name available
                            in inner tabular results file a tuple(TC, NC, GA).

        Return
        (dict)              Dictionary associating to each query name (with at
                            least one retrieved target name), the target row
                            which meets given GA threshold.

        Raise
        (KeyError)          In case no score has been set for a found query name
        """

        # Initailize hits dictionary associating queries to list of targets
        hits = dict()
        # Loop through each row in inner tabular results file
        for row in self.table:

            # Retrieve query name
            curr_qname = str(row['query_name'])
            # Retrieve bit score
            curr_bits = float(row['bit_score'])

            # Case no score is available for current hit (not significant)
            if curr_qname not in set(scores.keys()):
                # Go to next iteration
                continue

            # Retrieve scores for current query name
            tc, nc, ga = scores.get(curr_qname)

            # Case current bit score satisfies gathering threshold (GA)
            if curr_bits > ga:
                # Initialize set for current query name
                hits.setdefault(curr_qname, list())
                # Add current target name to hits set
                hits[curr_qname] += [row]

        # Return hits dictionary
        return hits


class Domtblout(Tblout):

    # Define number of columns
    width = 23

    # Define column names
    columns = {
        'target_name': 0,
        'query_name': 3,
        'e_value': 12,  # Use independent E-value
        'bit_score': 13,
        'alignment_beg': 17,
        'alignment_end': 18,
        'envelope_beg': 19,
        'envelope_end': 20
    }


# Unit testing
if __name__ == '__main__':

    # Define root directory
    ROOT_PATH = os.path.dirname(__file__) + '/../..'
    # Define temporary directory path
    TEMP_PATH = ROOT_PATH + '/tmp'
    # Example clusters directory
    EXAMPLE_PATH = TEMP_PATH + '/examples/MGYP*'

    # Define path to first cluster
    cluster_path = next(iglob(EXAMPLE_PATH))
    # Get first HMM model
    hmm_path = os.path.join(cluster_path, 'HMM')
    # Get FASTA file
    fasta_path = os.path.join(cluster_path, 'FA')

    # Define temporary output file
    out_path = NamedTemporaryFile(dir=os.getcwd(), delete=False, suffix='.sto').name
    # Log
    print('Making HMM alignment of {:s}'.format(out_path))

    # Define hmmalign script wrapper instance
    hmm_align = HMMAlign()

    # Try to run script
    try:
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Run hmmalign script
        hmm_align(hmm_path, fasta_path, out_path)
        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Log
        print('Generated HMM alignment in {:.0f}'.format(time_tot), end=' ')
        print('at {:s}'.format(out_path))

    except CalledProcessError as err:
        # Show error
        print('Error: hmmalign exited with code {}:'.format(err.returncode))
        print(err.stderr.strip())
        print(err.stdout.strip())
        # Raise error
        raise
