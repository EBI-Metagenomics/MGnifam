class Scheduler(object):
    """ Wrapper for Dask Cluster and Client

    Scheduler wraps Dask Cluster and Client, exposing methods for dinamically
    interact with remote/local clutser itself.
    """

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    # Update the cluster configuration, require workers
    def update(self, *args, min_jobs=0, max_jobs=100, **kwargs):
        raise NotImplementedError

    # Submit jobs, return results
    def submit(self, *args, **kwargs):
        # Sumbit computation, retrieve futures
        futures = self.client.submit(*args, **kwargs)
        # Retrieve results out of futures
        results = self.client.gather(futures)
        # Return results
        return results


class LocalScheduler(Scheduler):

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError


class LSFScheduler(Scheduler):

    # Constructor
    def __init__(self, *args, **kwargs):
        raise NotImplementedError
