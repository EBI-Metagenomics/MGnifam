# Dependencies
from src.hmm.hmmer import HMMER
from src.dataset import Fasta
from glob import iglob
import subprocess
import tempfile
import os
import re


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

    @classmethod
    def get_memory(cls, model_len, longest_seq=4e04, num_bytes=28):
        return model_len * longest_seq * num_bytes

    @classmethod
    def get_cpus(cls, max_memory, model_len, max_cpus=1, longest_seq=4e04, num_bytes=48):
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
        req_memory = cls.get_memory(model_len, longest_seq, num_bytes)
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


class Hits(object):
    """ Iterate through hmmsearch output

    This class parses a tabular output file from HMMSearch and retruns its
    attributes iteratively. However, not all the attributes are returned, only
    the ones needed in further computations.
    """
    def __init__(self, path):
        """ Constructor

        Args
        path (str)      Path to file storing tabular output results of hmmsearch
        """
        # Store input path
        self.path = path

    def open(self):
        # Open inner file
        self.file = open(self.path, 'r')
        # Return file reference
        return self.file

    def __enter__(self):
        # Open inner file
        self.open()
        # Return reference to self
        return self

    def close(self):
        # Close inner file
        self.file.close()
        # Remove inner file reference
        del self.file

    def __exit__(self, type, value, traceback):
        # Close and delete inner file buffer
        self.close()
        # Return itself
        return self

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
        (generator)                 List of rows: each row is a list, attributes
                                    in row list are expressed as
                                    tuple(non unique key, value)
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

    def get_iterator(self):
        """ Return iterator through tabular results file rows

        Return
        (generator)             Generator for loading file rows on the fly, as
                                returned by iterator static method

        Raise
        (FileNotFoundError)
        """
        # Ensure that inner file buffer exists
        if not hasattr(self, 'file'):
            # If not: raise new error
            raise FileNotFoundError(message=(
                'No file opened: inner file must be opened with either '
                '`open(...)` method or `with` statement!'
            ))
        # Return iterator over input tabular results file
        return self.iterator(self.file)

    def __iter__(self):
        return self.get_iterator()


# Implement sequence hits table parser
class SequenceHits(Hits):

    # Static columns dict(name: index)
    columns = {
        'target_name': 0,
        'query_name': 2,
        'e_value': 4,
        'bit_score': 5
    }

    def get_iterator(self):
        # Loop through each row in tabular file
        for row in super().get_iterator():
            # Retrieve row dictionary (each row value is a tuple (key, value))
            row = {key: row[index][1] for key, index in self.columns.items()}
            # Yield current row
            yield row

    # Compute score (TC, NC, GA) for each query
    def get_scores(self, e_value=0.001, min_bits=25.0, max_bits=999999.99):
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
        """

        # Initialize scores dictionary mapping query name to tuple(TC, NC, GA)
        scores = dict()

        # Loop through each row in current file
        for row in self:

            # Retrieve query name
            curr_qname = str(row['query_name'])
            # Retrieve e_value from current row
            curr_eval = float(row['e_value'])
            # Retrieve bit score from current row
            curr_bits = float(row['bit_score'])

            # Ensure that entry for current query name is set
            # TC is set extremely high (should include every score)
            # NC is set to lowest possible bit score
            # GA will be recomputed at the first iteration
            scores.setdefault(curr_qname, (max_bits, min_bits, min_bits))
            # Define current TC, NC, GA
            tc, nc, ga = scores[curr_qname]

            # Check if current e-value is below threshold one
            if curr_eval < e_value:
                # Check if current bit score is lower than TC
                if curr_bits < tc:
                    # Update tc, ensure it is greater than minimum bit score
                    tc = max(min_bits, curr_bits)

            # Otherwise: current e-value is greater or equal than threshold one
            else:
                # Check if current bit score is higher than threshold one
                if curr_bits > nc:
                    # Ensure bit score it is lower than TC
                    curr_bits = min(tc, curr_bits)
                    # Update nc, ensure it is greater than minimum bit score
                    nc = max(curr_bits, min_bits)

            # Update gathering threshold (GA): compute midpoint between NC and TC
            ga = nc + ((tc - nc) / 2)

            # Update scores for current query name
            scores[curr_qname] = (tc, nc, ga)

        # Return scores dict
        return scores

    @staticmethod
    def merge_scores(scores_list):
        """ Merge different score dictionaries

        Loop through each score dictionary in input list, merges TC, NC and GA
        for the retrieved query names.

        Args
        scores_list (list)          List of scores dictionaries

        Return
        (dict)                      Merged scores dictionary
        """
        # Initialize scores output dictionary
        scores = dict()

        # Loop through each score dictionary in input list
        for curr_scores in scores_list:
            # Loop through every query name in current scores
            for query_name in curr_scores.keys():
                # Retrieve current TC, NC, GA
                curr_tc, curr_nc, curr_ga = curr_scores.get(query_name)

                # Case this query name is not present in current scores dict
                if query_name not in scores.keys():
                    # Insert current scores
                    scores.setdefault(query_name, (curr_tc, curr_nc, curr_ga))
                    # Go to next iteration
                    continue

                # Retrieve previous TC NC GA
                prev_tc, prev_nc, prev_ga = scores.get(query_name)
                # Set new TC: minimum between current and previous
                curr_tc = min(curr_tc, prev_tc)
                # Set new NC: maximum between current and previous
                curr_nc = max(curr_nc, prev_nc)
                # Ensure that current NC is not higher than current TC
                curr_nc = min(curr_tc, curr_nc)
                # COmpute new gathering threshold (GA)
                curr_ga = curr_nc + ((curr_tc - curr_nc) / 2)

                # Store new triple (TC, NC, GA)
                scores[query_name] = (curr_tc, curr_nc, curr_ga)

        # Return scores output dictionary
        return scores

    # Retrieve entries according to scores (TC, NC, GA) for each query
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
                            least one retrieved target name), the set of target
                            names which meet GA threshold in given scores.

        Raise
        (KeyError)          In case no score has been set for a found query name
        """

        # Initailize hits dictionary associating queries to set of targets
        hits = dict()

        # Loop through each row in inner tabular results file
        for row in self:

            # Retrieve required fields form current row
            curr_qname = row['query_name']
            curr_tname = row['target_name']
            curr_eval = row['e_value']
            curr_bits = row['bit_score']

            # Retrieve scores for current query name
            tc, nc, ga = scores[curr_qname]

            # Case current bit score satisfies gathering threshold (GA)
            if curr_bits >= ga:
                # Initialize set for current query name
                hits.setdefault(curr_qname, set())
                # Add current target name to hits set
                hits.get(curr_qname).add(curr_tname)

        # Return hits dictionary
        return hits

    @staticmethod
    def merge_hits(hits_list):
        # Initialize merged hits dictionary
        hits = dict()

        # Loop through each hits dictionary in input list
        for hits in hits_list:
            # Loop through each key in current hits dictionary
            for query_name, target_names in hits.items():
                # Ensure there is a set for current query name
                hits.setdefault(query_name, set())
                # Add current target names
                hits.get(query_name).union(target_names)

        # Return merged hits dictionary
        return hits


# Implement domain hits table parser
class DomainHits(SequenceHits):

    # Static columns dict(name: index)
    columns = {
        'target_name': 0,
        'query_name': 3,
        'e_value': 12,
        'bit_score': 13
    }


# unit testing
if __name__ == '__main__':

    # Define path to root
    ROOT_PATH = os.path.join(os.path.dirname(__file__), '..', '..')
    # Define path to example clusters
    EXAMPLES_PATH = os.path.join(ROOT_PATH, 'tmp', 'examples', 'MGYP*')
    # Define path to UniProt dataset
    UNIPROT_PATH = os.path.join(ROOT_PATH, 'data_', 'uniprot', 'chunk*.fa.gz')
    # Define path to MGnifam dataset
    MGNIFAM_PATH = os.path.join(ROOT_PATH, 'data_', 'mgnify', 'chunk*.fa.gz')

    # Define new UniProt dataset
    uniprot = Fasta.from_str(UNIPROT_PATH)

    # Get first chunk in UniProt dataset
    chunk = uniprot[0]
    # Get longest sequence in chunk
    longest_seq, longest_len, chunk_size = chunk.get_longest()
    # Print longest sequence and chunk size
    print('Longest sequence found is of size {:d},'.format(longest_len), end=' ')
    print('found among other {:d} sequences'.format(chunk_size))
    print()

    # Retrieve first example cluster
    cluster_path = next(iglob(EXAMPLES_PATH))
    # Get chosen cluster name
    cluster_name = os.path.basename(cluster_path)
    # Show which cluster has been chosen
    print('Chosen cluster is {:s} at {:s}'.format(cluster_name, cluster_path))
    print()

    # Initailize tblout and domtblout files
    tblout_path = tempfile.NamedTemporaryFile(delete=False, suffix='.tblout').name
    domtblout_path = tempfile.NamedTemporaryFile(delete=False, suffix='.domtblout').name

    # Retrieve HMM from given cluster
    hmm_path = os.path.join(cluster_path, 'HMM')
    # Retrieve first chunk of UniProt
    chunk_path = chunk.path
    # Initialize new hmmsearch script
    hmm_search = HMMSearch()
    # Search HMM against first UniProt chunk
    hmm_search.run(
        hmm_path=hmm_path,
        target_path=chunk_path,
        tblout_path=tblout_path,
        domtblout_path=domtblout_path,
        seq_e=1e03,
        seq_z=chunk_size,
        num_cpus=1
    )

    # Print output paths
    print('Files retrieved from HMMSEARCH:')
    print('  {:s}'.format(tblout_path))
    print('  {:s}'.format(domtblout_path))
    print()

    # Test iterator for tblout (sequences)
    print('Retrieved tblout rows (sequence hits):')
    # Open tblout file
    with open(tblout_path, 'r') as tblout_file:
        # Define tblout iterator
        tblout_iter = Hits.iterator(tblout_file)
        # Loop through each row
        for i, row in enumerate(tblout_iter):
            # Print first three rows only
            if i > 3:
                continue
            # Print row
            print(row)
        # Print number of non printed hits
        print('...and {:d} other'.format(i - 3))
        print()

    # Test iterator for domtblout (domains)
    print('Retrieved domtblout rows (domain hits):')
    # Open domtblout file
    with open(domtblout_path, 'r') as domtblout_file:
        # Define tblout iterator
        domtblout_iter = Hits.iterator(domtblout_file)
        # Loop through each row
        for i, row in enumerate(domtblout_iter):
            # Print first three rows only
            if i > 3:
                continue
            # Print row
            print(row)
        # Print number of non printed hits
        print('...and {:d} other'.format(i - 3))
        print()

    # Test TBLOUT iterator (sequence hits)
    print('Retrieved tblout (sequences) hits:')
    # Read TBLOUT file
    with SequenceHits(tblout_path) as hits:
        # Loop through each hit
        for i, hit in enumerate(hits):
            # Print only the first 3 hits
            if i > 3:
                continue
            # Print current hit
            print(hit)
        # Print number of non printed hits
        print('...and {:d} other'.format(i - 3))
        print()

    # Test DOMTBLOUT iterator (domain hits)
    print('Retrieved domtblout (domains) hits:')
    # Read DOMTBLOUT file
    with DomainHits(domtblout_path) as hits:
        # Loop through each hit
        for i, hit in enumerate(hits):
            # Print only the first 3 hits
            if i > 3:
                continue
            # Print current hit
            print(hit)
        # Print number of non printed hits
        print('...and {:d} other'.format(i - 3))
        print()

    # Test sequneces scores retrieval
    with SequenceHits(tblout_path) as hits:
        # Retrieve scores
        scores = hits.get_scores(e_value=0.1)
        # Retrieve sequences
        sequences = hits.get_hits(scores)

        # Show sequence scores
        print('Retrieved sequence scores (e-value set to 0.1): ')
        # Loop through each query name
        for curr_qname, curr_scores in scores.items():
            # Print eithr query name and its scores
            print(curr_qname, curr_scores)

        # Show retrieved sequences
        print('Retrieved sequence hits: ')
        # Loop through each retrieved sequence
        for curr_qname, curr_tname in sequences:
            # Print either current query name and target name
            print(curr_qname, curr_tname)

    # Test sequneces scores retrieval
    with DomainHits(domtblout_path) as hits:
        # Retrieve scores
        scores = hits.get_scores(e_value=0.1)
        # Retrieve sequences
        domains = hits.get_hits(scores)

        # Show sequence scores
        print('Retrieved domain scores (e-value set to 0.1): ')
        # Loop through each query name
        for curr_qname, curr_scores in scores.items():
            # Print eithr query name and its scores
            print(curr_qname, curr_scores)

        # Show retrieved sequences
        print('Retrieved domain hits: ')
        # Loop through each retrieved sequence
        for curr_qname, curr_tname in domains:
            # Print either current query name and target name
            print(curr_qname, curr_tname)

    # Remove temporary files
    os.remove(tblout_path)
    os.remove(domtblout_path)
