"""HMM scan and HMM build

This library wraps HMMER functions for either building an Hidden Markov Model
(HMM) starting from a Multiple Sequence Alignment (MSA) and for searching a
HMM against a given dataset (e.g. UniProt).
"""


# Dependecies
import subprocess
import tempfile
import sys
import os


# Custom dependencies
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
    def run(self, msa, out_path, alphabet='amino'):

        # Initialize command line
        cmd = [*self.cmd]

        # Case alphabet type is not valid
        if alphabet not in self.alphabet:
            # Raise new error
            raise ValueError('Worng value for alphabet: {}'.format(alphabet))
        # Set type of alphabet used
        cmd += ['--{:s}'.format(alphabet)]

        # Set output file path
        cmd += [out_path]

        # Define new input MSA temporary file
        aln_file = tempfile.NamedTemporaryFile(suffix='.aln', delete=False)
        # Write out input MSA
        msa.to_aln(aln_file.name)
        # Set input MSA file in script
        cmd += [aln_file.name]

        # Run HMM build
        ran = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )

        # DEBUG Print executed command
        print('cmd:', ran)

        # Remove input file
        os.remove(aln_file.name)


class HMMSearch(HMMER):

    # Constructor
    def __init__(
        self,
        cmd=['/nfs/production/xfam/pfam/software/bin/hmmsearch'],
        env=os.environ.copy()
    ):
        # Call parent constructor
        super(HMMSearch, self).__init__(cmd=cmd, env=env)


# Unit testing
if __name__ == '__main__':

    # Define project root path
    ROOT_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
    # Path to test seed file
    SEED_PATH = os.path.join(ROOT_PATH, 'tmp', 'MGYP001224746368', 'SEED_trimmed')

    # Grab an example MSA
    msa = MSA().from_aln(SEED_PATH)
    # Debug
    print('Input MSA has shape: {}'.format(msa.aln.shape))

    # Define a new instance of hmm build script
    hmm_build = HMMBuild()
    # Run HMM build
    hmm_build.run(msa=msa, out_path=os.path.join(ROOT_PATH, 'tmp', 'hmm.out'))
