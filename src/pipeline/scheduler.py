class Scheduler(object):

    def __init__(self, cluster_type, **kwargs):
        # Save cluster type
        self.cluster_type = cluster_type
        # Save cluster arguments (as dictionary)
        self.cluster_kwargs = kwargs

    # Retrieve new cluster
    def get_cluster(self, **kwargs):
        # Override saved args with given kwargs (not persistently tough)
        kwargs = {**self.cluster_kwargs, **kwargs}
        # Instantiate cluster with given kwargs
        cluster = self.cluster_type(**kwargs)
        # Return the new cluster
        return cluster

    # Retrieve new client
    def get_client(self, **kwargs):
        # Just create new cluster
        cluster = self.get_cluster(**kwargs)
        # Generate new client over the created cluster
        client = Client(cluster)
        # return either the cluster and the client
        return cluster, client

    # Open
    def __enter__(self, **kwargs):
        # Return new cluster and new client
        return cluster, client

    # Close
    def __exit__(self):
        # TODO Close cluster
        # TODO Close client
        raise NotImplementedError
