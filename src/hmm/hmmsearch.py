# Dependencies
from src.hmm.hmmer import HMMER
from src.dataset import Fasta
from glob import iglob
import subprocess
import tempfile
import os


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
    def get_memory(cls, model_len, longest_seq=4e04, num_bytes=48):
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
        req_memory = cls.get_memory(model_len=model_len, longest_seq=longest_seq, num_bytes=num_bytes)
        # Loop from maximum number of CPUS to minimum (1)
        for num_cpus in range(max_cpus, 0, -1):
            # Check memory allocated to a single cpu
            cpu_memory = max_memory / num_cpus
            # Case allocated memory does not fit required
            if cpu_memory < req_memory:
                continue  # Skip iteration
            # Otherwise, return current number of cpus
            return num_cpus
        # If flow reaches this point, no cpu can be allocated
        return 0


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

    # Define model length
    model_len = 1000
    # Define longest sequence length
    longest_len = 4e04
    # Define maximum allowed memory per job
    max_memory = 8e09
    # Compute required memory per CPU
    req_memory = HMMSearch.get_memory(model_len, longest_len)
    # Compute maximum number of cores
    print(' '.join([
        'Memory required for a model of length {:d}'.format(model_len),
        'with a maximum sequence length of {:d}'.format(int(longest_len)),
        'is {:.2f} GB\n'.format(req_memory / 1e09)
    ]))

    # Compute max number of CPU with allowed memory 8 GB
    num_cpus = HMMSearch.get_cpus(max_memory=max_memory, model_len=model_len,
                                  max_cpus=8, longest_seq=longest_len)
    # Show retrieved number of CPUs
    print(' '.join([
        'Maximum number of CPUs allowed with maximum',
        'memory per job {:.2f} GB'.format(max_memory / 1e09),
        'is {:d}\n'.format(num_cpus)
    ]))

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

    # Define new MGnifam dataset
    mgnifam = Fasta.from_str(MGNIFAM_PATH)
    # Define MGnifam dataset chunk
    chunk = mgnifam[0]

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
        seq_z=900000,
        num_cpus=1
    )

    # Print output paths
    print('Files retrieved from HMMSEARCH:')
    print('  {:s}'.format(tblout_path))
    print('  {:s}'.format(domtblout_path))
    print()

    # Open TBLOUT file
    tblout = Tblout(tblout_path)
    # Print TBLOUT
    print('Tblout:')
    for row in tblout.table:
        print(row)

    # Open DOMTBLOUT file
    domtblout = Domtblout(domtblout_path)
    # Print DOMTBLOUT
    print('Domtblout:')
    for row in domtblout.table:
        print(row)

    # Remove either TBLOUT and DOMTBLOUT files
    os.remove(tblout_path)
    os.remove(domtblout_path)

    # # Test iterator for tblout (sequences)
    # print('Retrieved tblout rows (sequence hits):')
    # # Open tblout file
    # with open(tblout_path, 'r') as tblout_file:
    #     # Define tblout iterator
    #     tblout_iter = Hits.iterator(tblout_file)
    #     # Loop through each row
    #     for i, row in enumerate(tblout_iter):
    #         # Print first three rows only
    #         if i > 3:
    #             continue
    #         # Print row
    #         print(row)
    #     # Print number of non printed hits
    #     print('...and {:d} other'.format(i - 3))
    #     print()
    #
    # # Test iterator for domtblout (domains)
    # print('Retrieved domtblout rows (domain hits):')
    # # Open domtblout file
    # with open(domtblout_path, 'r') as domtblout_file:
    #     # Define tblout iterator
    #     domtblout_iter = Hits.iterator(domtblout_file)
    #     # Loop through each row
    #     for i, row in enumerate(domtblout_iter):
    #         # Print first three rows only
    #         if i > 3:
    #             continue
    #         # Print row
    #         print(row)
    #     # Print number of non printed hits
    #     print('...and {:d} other'.format(i - 3))
    #     print()
    #
    # # Test TBLOUT iterator (sequence hits)
    # print('Retrieved tblout (sequences) hits:')
    # # Read TBLOUT file
    # with SequenceHits(tblout_path) as hits:
    #     # Loop through each hit
    #     for i, hit in enumerate(hits.iterate()):
    #         # Print only the first 3 hits
    #         if i > 3:
    #             continue
    #         # Print current hit
    #         print(hit)
    #     # Print number of non printed hits
    #     print('...and {:d} other'.format(i - 3))
    #     print()
    #
    # # Test DOMTBLOUT iterator (domain hits)
    # print('Retrieved domtblout (domains) hits:')
    # # Read DOMTBLOUT file
    # with DomainHits(domtblout_path) as hits:
    #     # Loop through each hit
    #     for i, hit in enumerate(hits.iterate()):
    #         # Print only the first 3 hits
    #         if i > 3:
    #             continue
    #         # Print current hit
    #         print(hit)
    #     # Print number of non printed hits
    #     print('...and {:d} other'.format(i - 3))
    #     print()
    #
    # # Test sequneces scores retrieval
    # with SequenceHits(tblout_path) as hits:
    #
    #     # Retrieve scores
    #     scores = hits.get_scores(e_value=0.1)
    #     # Show sequence scores
    #     print('Retrieved sequence scores (e-value set to 0.1): ')
    #     # Loop through each query name
    #     for curr_qname, curr_scores in scores.items():
    #         # Print eithr query name and its scores
    #         print(curr_qname, curr_scores)
    #
    #     # Retrieve sequences
    #     sequences = hits.get_hits(scores)
    #     # Show retrieved sequences
    #     print('Retrieved sequence hits: ')
    #     # Loop through each retrieved sequence
    #     for curr_qname, curr_tname in sequences:
    #         # Print either current query name and target name
    #         print(curr_qname, curr_tname)
    #
    # # Test sequneces scores retrieval
    # with DomainHits(domtblout_path) as hits:
    #
    #     # Retrieve scores
    #     scores = hits.get_scores(e_value=0.1)
    #     # Show sequence scores
    #     print('Retrieved domain scores (e-value set to 0.1): ')
    #     # Loop through each query name
    #     for curr_qname, curr_scores in scores.items():
    #         # Print eithr query name and its scores
    #         print(curr_qname, curr_scores)
    #
    #     # Retrieve sequences
    #     domains = hits.get_hits(scores)
    #     # Show retrieved sequences
    #     print('Retrieved domain hits: ')
    #     # Loop through each retrieved sequence
    #     for curr_qname, curr_tname in domains:
    #         # Print either current query name and target name
    #         print(curr_qname, curr_tname)
    #
    # # Remove temporary files
    # os.remove(tblout_path)
    # os.remove(domtblout_path)
    #
    # # Initialize scores dictionary
    # scores = dict()
    # # Define UniProt iterator
    # uniprot_iter = iglob(UNIPROT_PATH)
    # # Loop through first 3 example HMM models
    # for i, chunk_path in enumerate(uniprot_iter):
    #     # Initailize tblout and domtblout files
    #     tblout_path = tempfile.NamedTemporaryFile(delete=False).name
    #     domtblout_path = tempfile.NamedTemporaryFile(delete=False).name
    #
    #     # Make HMM search
    #     hmm_search.run(
    #         hmm_path=hmm_path,
    #         target_path=chunk_path,
    #         tblout_path=tblout_path,
    #         domtblout_path=domtblout_path,
    #         seq_e=1000,
    #         seq_z=1e06,
    #         num_cpus=1
    #     )
    #
    #     # Retrieve scores
    #     with SequenceHits(tblout_path) as hits:
    #         # Retrieve scores
    #         scores = hits.merge_scores([
    #             scores,  # Previously stored scores
    #             hits.get_scores(e_value=1)  # Newly retrieved scores
    #         ])
    #
    #     # Print scores
    #     print('Scores:')
    #     for query_name, score in scores.items():
    #         print(query_name, score)
    #     print()
    #
    #     # Remove temporary files
    #     os.remove(tblout_path)
    #     os.remove(domtblout_path)
    #
    #     # Early stopping condition
    #     if i >= 3:
    #         break
