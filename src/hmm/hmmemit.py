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
    def run(self, out_path='', num_samples=None, emit_type=None):
        # Initialize command
        cmd = self.cmd
        # Eventually add output file path
        cmd += ['-o', out_path] if out_path else []
        # Eventually add en emission type
        cmd += [emit_type] if emit_type else []
        # Ensure that each input is string
        cmd = [str(x) for x in cmd]

        # Run process, return result object
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )
