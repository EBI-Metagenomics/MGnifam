# Dependencies
import json
import sys
import os


class Log(object):

    # Cosntructor
    def __init__(self, log_dict=dict(), log_path=''):
        # Store attributes
        self.log_dict = log_dict
        self.log_path = log_path

    # Wrapper for update method
    def __call__(self, *args, **kwargs):
        # Call update method
        return self.update(*args, **kwargs)

    # Update method
    def update(self, log_dict=dict(), to_file=True):
        # Update inner log dictionary
        self.log_dict = {**self.log_dict, **log_dict}
        # Case log file path is set
        if self.log_path and to_file:
            # Dump new log to file (overwrites the previous one)
            self.to_file()

    # Log to file
    def to_file(self):
        # Open output file
        with open(self.log_path, 'w') as log_file:
            # Dump log dictionary to file
            json.dump(log_dict, log_file, indent=2)


class Pipeline(object):

    # (Abstract) Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    # Wrapper for run method
    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    # (Abstract) Run method
    def run(self, *args, **kwargs):
        raise NotImplementedError
