# Dependencies
import json
import sys
import os

# Dask
from dask.distributed import Client


class Log(object):

    # Cosntructor
    def __init__(self, log_dict=dict(), log_path=''):
        # Store attributes
        self.log_dict = log_dict
        self.log_path = log_path

    # Override squared brackets
    def __getitem__(self, key):
        # Define underlying dictionary keys
        keys = set(self.log_dict.keys())
        # Case key does not exist in underlying dictionary
        if key not in keys:
            # Raise key error
            raise KeyError('Key not found in underlying log dictionary')
        # Return key in log dictionary
        return self.log_dict[key]

    # Override squared brackets
    def __setitem__(self, key, item):
        # Set item in inner log dictionary
        self.log_dict[key] = item

    # Handle attribute not found: search in log dict
    def __getattr__(self, key):
        # Wrapper for __getitem__
        return self[key]

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
            json.dump(self.log_dict, log_file, indent=2)


class Pipeline(object):

    # Constructor
    def __init__(self, cluster_type, cluster_kwargs):
        # Store dask arguments
        self.cluster_type = cluster_type
        self.cluster_kwargs = cluster_kwargs

    # Make Dask cluster
    def get_cluster(self, cluster_kwargs=dict()):
        # Update default cluster kwargs
        cluster_kwargs = {**self.cluster_kwargs, **cluster_kwargs}
        # Make new cluster
        return self.cluster_type(**cluster_kwargs)

    # Make Dask client
    def get_client(self, *args, **kwargs):
        # Define new cluster
        cluster = self.get_cluster(*args, **kwargs)
        # Define new client
        client = Client(client)
        # Return client containing cluster
        return cluster, client

    # Close Dask cluster and client
    @staticmethod
    def close_client(cluster, client):
        # Close cluster first
        cluster.close()
        # Then, close client
        client.close()

    # Wrapper for run method
    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    # (Abstract) Run method
    def run(self, *args, **kwargs):
        raise NotImplementedError
