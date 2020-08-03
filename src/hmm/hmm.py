# Dependencies
from src.utils import open_file
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
    def get_name(cls, string):
        # Check if input string matches name string
        is_name = re.search(r'^NAME\s+(\S+)', string)
        # Case string is not name
        if not is_name:
            return None
        # Otherwise, retrieve and return name
        return is_name.group(1)

    @classmethod
    def get_length(cls, string):
        # Check if input string matches length string
        is_length = re.search(r'^LENG[\s]+(\d+)', string)
        # Case string is not length
        if not is_length:
            return None
        # Otherwise, return length string
        return int(is_length.group(1))

    @classmethod
    def get_alphabet(cls, string):
        # Check if input string matches alphabet string
        is_alphabet = re.search(r'^ALPH\s+(\S+)', string)
        # Case string is not alphabet
        if not is_alphabet:
            return None
        # Otherwise, return alphabet string
        return is_alphabet.group(1)

    @classmethod
    def from_file(cls, path):
        # Define new dict containing default HMM parameters
        params = {'name': '', 'length': 0, 'alphabet': 'amino'}
        # Open input file
        file = open(path, 'r')
        # Loop through each line in file
        for line in file:
            # Get name
            name = cls.get_name(line)
            # Get length
            length = cls.get_length(line)
            # Get alphabet
            alpha = cls.get_alphabet(line)
            # Update params fields
            params['name'] = params['name'] if name is None else name
            params['length'] = params['length'] if length is None else length
            params['alphabet'] = params['alphabet'] if alpha is None else alpha
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
