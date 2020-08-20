# Dependencies
from src.hmm.hmmer import HMMER
import subprocess
import os


class HMMEmit(HMMER):

    # Define emit type constants (check hmmemit )
    EMIT_ALIGN = '-a'
    EMIT_CONSENSUS = '-c'
    EMIT_FANCIER = '-C'
    EMIT_UNALIGN = '-p'

    # Constructor
    def __init__(self, cmd=['hmmemit'], env=os.environ):
        # Call parent constructor
        super().__init__(cmd=cmd, env=env)

    # Run
    def run(self, model_path, out_path='', num_samples=None, emit_type=None):
        """ Run hmmemit

        Args
        model_path (str)        Path to input HMM model file
        out_path (str)          Path to output file
        num_samples (int)       Number of samples to retrieve, defaults to 10
        emit_type (str)         Type of emission to perform

        Return
        CompletedProcess        Results retrieved successfully running script

        Raise
        (FileNotFoundError)     In case input file does not exist
        (CalledProcessError)    In case errors occurred while running script
        """
        # Initialize command
        cmd = self.cmd
        # Eventually add output file path
        cmd += ['-o', out_path] if out_path else []
        # Eventually add en emission type
        cmd += [emit_type] if emit_type else []
        # Add mandatory input model file
        cmd += [model_path]
        # Ensure that each input is string
        cmd = [str(c) for c in cmd]

        # Check input file
        if not os.path.isfile(model_path):
            # Raise exception
            FileNotFoundError(' '.join([
                'could not find input HMM model file',
                'at {:s}'.format(model_path)
            ]))

        # Run process, return result object
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )
