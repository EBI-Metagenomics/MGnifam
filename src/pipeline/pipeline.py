# Dependencies
import json
import sys
import os


class Pipeline(object):

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    # Run method
    def run(self, *args, **kwargs):
        raise NotImplementedError

    # Call method (wrap run method)
    def __call__(self, *args, **kwargs):
        # Call run method
        return self.run(*args, **kwargs)
