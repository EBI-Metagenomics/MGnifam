# Dependencies
from dask.distributed import Client
from dask_jobqueue import LSFCluster
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

# Main
if __name__ == '__main__':

    # Make new cluster
    cluster = LSFCluster(walltime='01:00', memory='1GB', cores=1, processes=1)
    # Make new client
    client = Client(cluster)

    # Define interval boundaries (count to 100)
    count_beg, count_end = 0, 100
    # Define batch size
    batch_size = 10

    # Define minimum and maximum number of jobs to use
    min_jobs, max_jobs = 0, 1

    # Request for 10 jobs
    # cluster.adapt(minimum=min_jobs, maximum=max_jobs)
    cluster.adapt()

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
        futures.append(client.submit(test, batch_beg, batch_end, delay=1))
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
