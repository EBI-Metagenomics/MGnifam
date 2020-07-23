class Scheduler(object):

    # Constructor
    def __init__(self, cluster_type, **kwargs):
        # Save cluster type
        self.cluster_type = cluster_type
        # Save cluster kwargs
        self.cluster_kwargs = kwargs
        # Define inner cluster
        self.cluster = self.cluster_type(**self.cluster_kwargs)
        # Define inner CLient
        self.client = Client(self.cluster)

    # Stop connection to cluster
    def close(self):
        # Close inner cluster
        self.cluster.close()
        # Close inner client
        self.client.close()

    # Open
    def __enter__(self, **kwargs):
        # Just return self
        return self

    # Close
    def __exit__(self):
        # Just call close method
        self.close()
