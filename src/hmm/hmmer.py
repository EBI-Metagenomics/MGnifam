# Dependencies
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
from glob import iglob
from time import time
import subprocess
import os


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
    def __init__(self, cmd=['hmmbuild'], env=os.environ.copy()):
        # Call parent constructor
        super(HMMBuild, self).__init__(cmd=cmd, env=env)

    # Run HMM build script
    def run(
        self, msa_path, hmm_path, name=None, log_path='/dev/null',
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
        cmd += [hmm_path]
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


class HMMAlign(HMMER):

    # Constructor
    def __init__(self, cmd=['hmmalign'], env=os.environ.copy()):
        # Call parent constructor
        super().__init__(cmd=cmd, env=env)

    # Run HMM align
    def run(self, hmm_path, fasta_path, out_path):
        """Align an HMM model with a FASTA dataset

        Args
        hmm_path (str)          Path to HMM model file
        fasta_path (str)        Path to FASTA target sequences file
        out_path (str)          Path where to store output alignment file

        Return
        (CompletedProcess)      Result object representing a successful run

        Raise
        (CalledProcessError)    Error related to command line input
        (Error)                 Unknown error
        """
        # Define command
        cmd = self.cmd
        # Set mandatory command line options
        cmd += ['-o', out_path]  # Output file path
        cmd += [hmm_path]  # HMM model file path
        cmd += [fasta_path]   # Target dataset file path

        # Execute script
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )


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
