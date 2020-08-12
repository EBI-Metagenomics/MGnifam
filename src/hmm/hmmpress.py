# Dependencies
from src.hmm.hmmer import HMMER
from glob import iglob
import subprocess
import argparse
import os


class HMMPress(HMMER):

    # Constructor
    def __init__(self, cmd=['hmmpress'], env=os.environ):
        # Call parent constructor
        super().__init__(cmd=cmd, env=env)

    # Run script
    def run(self, hmm_path, overwrite=False):
        """ Run hmmpress script

        Args
        hmm_path (str)          Path to HMM model/library
        overwrite (str)          Whether generated .h3[mifp] files can be
                                overwritten

        Return
        (CompletedProcess)      Script output

        Raise
        (CalledProcessError)    When script did not run properly,
                                for example if generated files can not be
                                overwritten and are already in place
        """
        # Intialize commands
        cmd = self.cmd
        # Eventually force files overwrite
        cmd += ['-f'] if overwrite else []
        # Add mandatory parameters: HMM path
        cmd += [hmm_path]

        # Run hmmpress, return CompletedProcess instance
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )


# Run hmmpress from command line
if __name__ == '__main__':

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Run hmmpress')
    # Initialize hmmscan command
    parser.add_argument(
        '--cmd', nargs='+', type=str, default=['hmmpress'],
        help='Path to hmmpress executable'
    )
    # Define HMM path
    parser.add_argument(
        '--hmm_path', type=str, required=True,
        help='Path to input HMM model/library'
    )
    # Set overwrite flag
    parser.add_argument(
        '--overwrite', type=bool, required=False,
        help='Whether generated files can be automatically overwritten'
    )

    # Retrieve arguments
    args = parser.parse_args()

    # Try running hmmscan
    try:
        # Check if file exists
        if not os.path.isfile(args.hmm_path):
            # Raise file not found error
            raise FileNotFoundError(' '.join([
                'Error: unable to find input HMM model/library file',
                'at {:s}'.format(args.hmm_path)
            ]))

        # Make new HMMPress instance
        hmmpress = HMMPress(cmd=args.cmd, env=os.environ)
        # Run HMMPress on given HMM model/library
        hmmpress(hmm_path=args.hmm_path, overwrite=args.overwrite)

        # Retrieve output directory
        hmm_dir = os.path.dirname(args.hmm_path)
        # Search for output files
        h3m_path = next(iglob(hmm_dir + '/*.h3m'), '')
        h3i_path = next(iglob(hmm_dir + '/*.h3i'), '')
        h3f_path = next(iglob(hmm_dir + '/*.h3f'), '')
        h3p_path = next(iglob(hmm_dir + '/*.h3p'), '')

        # Verbose
        print('Generated .h3[mifp] files:')
        print('- h3m_path: {:s}'.format(h3m_path))
        print('- h3i_path: {:s}'.format(h3i_path))
        print('- h3f_path: {:s}'.format(h3f_path))
        print('- h3p_path: {:s}'.format(h3p_path))

    except subprocess.CalledProcessError as err:
        # Raise new exception
        raise Exception('\n'.join([
            'Error: hmmpress returned code {}:'.format(err.returncode),
            err.stderr.strip()
        ]))
