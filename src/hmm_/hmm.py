# Dependencies
import sys
import os
import re


class HMM(object):

    # Constructor
    def __init__(self, name='', length=0, alphabet='amino'):
        # HMM attributes
        self.name = name
        self.length = length
        self.alphabet = alphabet

    @classmethod
    def from_file(cls, path):
        # Define new dict containing default HMM parameters
        params = {'name': '', 'length': 0, 'alphabet': 'amino'}
        # Open input file
        file = open_file(path, 'r')
        # Loop through each line in file
        for line in file:
            # Check if line is HMM name
            is_name = re.search(r'^NAME\s+(\S+)', line)
            # Case line is HMM name
            if is_name:
                # Store name
                params['name'] = str(is_name.group(1))
                # Go to next iteration
                continue
            # Check if line is HMM length
            is_length = re.search(r'^LENG\s+(\d+)', line)
            # Case line is HMM length
            if is_length:
                # Store HMM length
                params['length'] = int(is_length.group(1))
                # Go to next iteration
                continue
            # Check if line is HMM alphabet
            is_alphabet = re.search(r'^ALPH\s+(\S+)', line)
            # Case line is HMM alphabet
            if is_alphabet:
                # Store HMM alphabet
                params['alphabet'] = str(is_alphabet.group(1))
                # Go to next iteration
                continue
        # Close file buffer
        file.close()
        # Return new HMM instance with retrieved params
        return cls(**params)

    @classmethod
    def concat(cls, in_paths, out_path):
        """Create HMM library
        Create an HMM library by concatenating multiple HMM files

        Args
        in_paths (list)     List of HMM files/libraries to concatenate
        out_path (str)      Path to resulting HMM library file

        Raise
        (FileNotFoundError) Error when trying to open one input file or when
                            creating the output file
        """
        # Open output file
        with open(out_path, 'w') as out_file:
            # Loop through each input path
            for in_path in in_paths:
                # Open input path
                in_file = open(in_path, 'r')
                # Loop through every line in input file
                for line in in_file:
                    # Write out line in output file
                    out_file.write(line)
                # Close input file
                in_file.close()
