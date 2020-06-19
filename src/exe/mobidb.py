"""Mobidb wrapper
Handles Mobidb executrable script
"""

# Dependencies
import subprocess
import logging as log
import os


class MobiDB(object):

    # Executable script (static)
    args = ['python3', 'mobidb_lite.py', '-o', '{out_path:s}', '{in_path:s}']

    @classmethod
    def run(cls, run_in, run_out='/dev/stdout', cwd='./', env=os.environ.copy()):
        # Debug
        log.debug('MobiDB lite input path: ' + str(run_in))
        log.debug('MobiDB lite output path: ' + str(run_out))
        # Get arguments
        args = cls.args.copy()
        # Set input and output paths
        args[3] = args[3].format(run_out)
        args[5] = args[5].format(run_in)
        # Run MobiDB Lite through a subprocess
        out = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # Set output encoding
            env=env,  # Set environment variables
            cwd=cwd,  # Set working directory
            args=args  # Set script arguments
        )
        # Debug output
        log.debug('MobiDB lite console output: ' + str(out))
        # Return output
        return out
        raise NotImplementedError
