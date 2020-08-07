# Dependencies
import dask.distributed as dd
import dask_jobqueue as dj
from time import time, sleep

# Define test function
def test(beg, end, delay=1):
    # Verbose
    print('This is a test function: ', end=' ')
    print('counts from {:d} to {:}', end=' ')
    print('with a delay of {:d} seconds')
    # Initialize return array
    numbers = list()
    # Loop through every integer between begin and end
    for i in range(beg, end):
        # Store current index
        numbers.append(i)
        # Wait for given delay
        sleep(delay)
    # Return list of numbers
    return numbers


# Just a wrapper for a LSFCluster with scalable cores
class LSFCluster(dj.LSFCluster):

    # Constructor
    def __init__(self, memory, cores, processes, walltime='00:30', **kwargs):
        # Call parent contructort
        super().__init__(memory=memory, cores=cores, processes=processes, walltime=walltime)
        # Save number of cores
        self.cores = cores
        # DEBUG
        print('Number of cores:', self.cores)
        print('Job script:', self.job_script())

    # Adapt: require this amount of cores per job
    def adapt(self, minimum, maximum, cores=1):
        raise NotImplementedError


# Main
if __name__ == '__main__':

    # Make new cluster
    # cluster = dj.LSFCluster(walltime='00:30', memory='1GB', cores=1, processes=1)
    cluster = LSFCluster(walltime='00:30', memory='1GB', cores=2, processes=1)
    # Make new client
    client = dd.Client(cluster)

    # Define interval boundaries (count to 100)
    count_beg, count_end = 0, 100
    # Define batch size
    batch_size = 10

    # Define minimum and maximum number of jobs to use
    min_jobs, max_jobs = 0, 1

    # Request for 10 jobs
    # cluster.adapt(minimum=min_jobs, maximum=max_jobs)
    cluster.scale(n=10, jobs=10, memory='4GB', cores=5)

    # Define start time
    time_beg, time_end = time(), 0.0
    # Initialize futures container
    futures = list()
    # Loop through each batch start index
    for i in range(count_beg, count_end, batch_size):
        # Define batch start index
        batch_beg = i
        # Define batch end index
        batch_end = i + min(batch_size, count_end)
        # Make new future
        futures.append(client.submit(test, batch_beg, batch_end, delay=3))
    # Retrieve results
    results = client.gather(futures)
    # Update end time
    time_end = time()
    time_tot = time_end - time_beg

    # Summary
    print('Took {:.0f} seconds to count'.format(time_tot), end=' ')
    print('from {:d} to {:d}'.format(count_beg, count_end))
    # Go through each result array
    for i in range(len(results)):
        # Print current result
        print('{:03d}-th result:'.format(i+1), end=' ')
        print(' '.join([str(j) for j in results[i]]))

    # Close client
    client.close()
    # Close cluster
    cluster.close()
