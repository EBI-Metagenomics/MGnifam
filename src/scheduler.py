# Dependencies
from dask.distributed import Client, LocalCluster
from distributed.utils import parse_bytes
from dask_jobqueue import LSFCluster
from time import time, sleep
import dask

# Setup Dask
dask.config.set({"distributed.comm.timeouts.connect": "30s"})


# Distributed/local jobs handler
class Scheduler(object):

    # Abstract constructor
    def __init__(self, min_cores=1, max_cores=1, min_memory='1 GB', max_memory='1 GB'):
        # Store bounaries
        self.min_cores = min_cores
        self.max_cores = max_cores
        self.min_memory = min_memory
        self.max_memory = max_memory
        # Initialize both client and cluster
        self._client = None
        self._cluster = None

    @property
    def cluster(self):
        return self._cluster

    @property
    def client(self):
        return self._client

    @property
    def max_cores(self):
        return self._max_cores

    @max_cores.setter
    def max_cores(self, max_cores):
        self._max_cores = max_cores

    @property
    def min_cores(self):
        return self._min_cores

    @min_cores.setter
    def min_cores(self, min_cores):
        self._min_cores = min_cores

    @property
    def max_memory(self):
        return self._max_memory

    @max_memory.setter
    def max_memory(self, max_memory):
        self._max_memory = max_memory

    @property
    def min_memory(self):
        return self._min_memory

    @min_memory.setter
    def min_memory(self, min_memory):
        self._min_memory = min_memory

    # Require n jobs with given cores and memory characteristics
    def adapt(self, minimum, maximum, cores=1, memory='1 GB', **kwargs):

        # Check if given memory is greater than maximum allowed
        if parse_bytes(memory) > parse_bytes(self.max_memory):
            raise MemoryError(' '.join([
                'could not allocate {:s} of memory,'.format(memory),
                'maximum allowed is {:s}'.format(self.max_memory)
            ]))

        # Check if given memory is lower than minimum allowed
        if parse_bytes(memory) < parse_bytes(self.min_memory):
            raise MemoryError(' '.join([
                'could not allocate {:s} of memory,'.format(memory),
                'minimum allowed is {:s}'.format(self.min_memory)
            ]))

        # Check if number of cores is greater than maximum allowed
        if cores > self.max_cores:
            raise Exception(' '.join([
                'could not allocate {:d} cores,'.format(cores),
                'maximum allowed is {:d}'.format(self.max_cores)
            ]))

        # Check if number of cores is lower than minimum allowed
        if cores < self.min_cores:
            raise Exception(' '.join([
                'could not allocate {:d} cores,'.format(cores),
                'minimum allowed is {:d}'.format(self.min_cores)
            ]))

    # Close either cluster and client
    def close(self):
        # Stop client
        self.client.close()
        # Stop cluster
        self.cluster.close()
        # Remove both cluster and client
        self._cluster, self._client = None, None


# LSF distributed jobs handler
class LSFScheduler(Scheduler):

    # Constructor
    def __init__(self, min_cores=1, max_cores=1, min_memory='1 GB', max_memory='1 GB', processes=1, walltime='02:00', **kwargs):
        # Call parent constructor
        super().__init__(min_cores=min_cores, max_cores=max_cores, min_memory=min_memory, max_memory=max_memory)
        # Define cluster default parameters
        self.cluster_kwargs = {**{
            'memory': max_memory,
            'cores': min_cores,
            'processes': processes,
            'walltime': walltime
        }, **kwargs}

    # Define adapt method
    def adapt(self, minimum, maximum, **kwargs):
        # Merge kwargs with default kwargs
        kwargs = {**self.cluster_kwargs, **kwargs}
        # Call parent adapt method (check values)
        super().adapt(minimum, maximum, **kwargs)
        # Make new cluster
        self._cluster = LSFCluster(**kwargs)
        # Make client
        self._client = Client(self._cluster)
        # Adapt cluster
        self._cluster.adapt(minimum=minimum, maximum=maximum)
        # Return client reference
        return self.client


# Local jobs handler
class LocalScheduler(Scheduler):

    # Constructor
    def __init__(self, min_cores=1, max_cores=1, min_memory='1 GB', max_memory='1 GB', processes=1, **kwargs):
        # Call parent constructor
        super().__init__(min_cores=min_cores, max_cores=max_cores, min_memory=min_memory, max_memory=max_memory)
        # Define cluster default parameters
        self.cluster_kwargs = {**{
            'threads_per_worker': min_cores,
            'processes': processes
        }, **kwargs}

    def adapt(self, minimum, maximum, **kwargs):
        # Merge kwargs with default kwargs
        kwargs = {**self.cluster_kwargs, **kwargs}
        # Call parent adapt method (check values)
        super().adapt(minimum, maximum, **kwargs)
        # Make new cluster
        self._cluster = LocalCluster(**kwargs)
        # Make client
        self._client = Client(self._cluster)
        # Adapt cluster
        self._cluster.adapt(minimum=minimum, maximum=maximum)
        # Return client reference
        return self.client


# Unit testing
if __name__ == '__main__':

    # Define test function
    def count_beg_end(beg, end, delay=1):
        # Verbose
        print('This is a test function: ', end=' ')
        print('it counts from {:d} to {:d}'.format(beg, end), end=' ')
        print('with a delay of {:d} seconds')
        # Loop through every integer between begin and end
        for i in range(beg, end, 1):
            # Wait for given delay
            sleep(delay)

    # Define scheduler
    scheduler = LSFScheduler(
        # Set cores boudaries
        min_cores=1,
        max_cores=3,
        # Set memory boundaries
        min_memory='1 GB',
        max_memory='3 GB',
        # Define processes (jobs) per worker
        processes=1
    )
    # scheduler = LocalScheduler(min_cores=1, max_cores=3, min_memory='1 GB', max_memory='3GB')

    # Initialize timers
    time_beg, time_end = time(), 0.0
    # Adapt
    with scheduler.adapt(10, 10, cores=1, memory='1 GB') as client:
        # Define futures
        futures = list()
        # Do ten (10) loops
        for i in range(0, 100, 10):
            # Count to 100
            future = client.submit(count_beg_end, beg=i, end=i+10, delay=2)
            # Store future
            futures.append(future)

        # Show number of futures
        print('Number of futures to wait: {:d}'.format(len(futures)))

        # Gather results
        client.gather(futures)
        # Update timers
        time_end = time()
        # Verbose
        print('It took {:.0f} seconds to count to 100 ten times'.format(time_end - time_beg))

        # Verbose
        print('Closing connection...', end='')
        # Wait for half a minute
        sleep(30)

    # Show connection closed
    print('done!')
    # Wait for another half a minute
    sleep(30)
