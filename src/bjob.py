"""
Handle job submitting to LSF
"""


# Dependencies
import subprocess
import time
import re


class Bjob(object):

    def __init__(self, id=None, status=None, out_path='/dev/null', err_path='/dev/null'):
        self.id = id
        self.curr_status = status
        self.out_path = out_path
        self.err_path = err_path

    def get_status(self):
        """Retrieve job status
        Check status onnly if it is RUN or PEND, otherwise, returns the already
        set status

        Return
        (str)   Job status, which can be 'RUN', 'PEND', 'DONE', 'EXIT'
        """
        # If not running, return current status
        if self.curr_status not in set(['RUN', 'PEND']):
            return self.curr_status
        # Otherwise, retrieve status using bjobs
        self.curr_status = self.status(self.id)
        # Return retrieved status
        return self.curr_status

    def is_running(self):
        """Define if current job is running

        Return
        (bool)  Whether job is running or not
        """
        # Get current job status
        status = self.get_status()
        # Return bool if status is running or pending
        return (status in set(['RUN', 'PEND']))

    def is_done(self):
        """States wether current job is done

        Return
        (bool)  True if job is DONE, false otherwise
        """
        return self.get_status() == 'DONE'

    @classmethod
    def run(cls, args, out_path='/dev/null', err_path='/dev/null', queue=None, walltime=None, cpu=1, memory='2 GB', verbose=False):
        """Run a job
        Run a job by submitting it to LSF scheduler: this method takes as input
        subprocess arguments and feeds them to bsub command, along with other
        optional arguments

        Args
        args (list):            List of subrpocess arguments
        out_path (str)          Path to output log file
        err_path (str)          Path to error log file
        queue (str)             On which queue command must be run
        walltime (str)          HH:MM a job can stay alive
        cpu (int)               Number of CPUs needed for that job
        memory (str)            How much memory is needed for that job
        verbose (bool)          Whether to print verbose log or not

        Return
        (Bjob):                 Instance of generated bjob

        Raise
        (CalledProcessError)    Unable to run command
        """
        # Verbose
        if verbose:
            print('Starting new job...')

        # Reformat code to make it run on LSF
        wrap = ['bsub', '-o', out_path, '-e', err_path]
        # TODO Set queue
        # TODO Set walltime
        # TODO Set number of CPU
        # TODO Set memory

        # Add user define arguments
        wrap += args

        # Run command on LSF
        ran = subprocess.run(
            check=True,  # Eventually raise error
            capture_output=True,  # Retain output
            encoding='utf-8',  # Encoding
            args=args  # bsub arguments
        )

        # Verbose
        if verbose:
            print('Started `{:s}`'.format(str(ran.args)))
            print('with stdout {:s}'.format(ran.stdout))
            print('     stderr {:s}'.format(ran.stderr))

        # Retrieve job id
        job_id = cls.id_from_string(ran.stdout)
        # Retrieved job status, running by default
        job_status = 'RUN'
        # Case id is not defined
        if job_id is None:
            # Turn job status to exit
            job_status='EXIT'

        # Return retrieved Bjob instance with give job id
        return cls(id=job_id, status=job_status, out_path=out_path, err_path=err_path)

    @classmethod
    def check(cls, bjobs, delay=30, verbose=False):
        """Check list of jobs
        Checks the given list of jobs (or a single one) until all are finished

        Args
        bjobs (list)        List of bjobs whose status must be checked
        delay (int)         Number of seconds between a check and the other
        verbose (bool)      Wether to return verbose output or not
        """
        # Start looping
        while True:

            # Initialize counter of running jobs
            num_running = 0
            # Loop through input list of jobs
            for bjob in bjobs:
                # Case job is not running
                if not bjob.is_running():
                    continue  # Go to next iteration
                # Signal that job is still running
                num_running += 1

            # Verbose
            if verbose:
                print('There are {:d} jobs'.format(num_running), end=' ')
                print('which are still running')

            # Check if at least one process is still running
            if num_running < 1:
                # Verbose
                if verbose:
                    print('No more running jobs')
                # Exit method
                return

            # Wait for given delay
            time.sleep(delay)

    @classmethod
    def status(cls, job_id, verbose=False):
        # Retrieve command output
        ran = subprocess.run(
            check=True,
            capture_output=True,
            encoding='utf-8',
            args=['bjobs', '-noheader', '-a', job_id]
        )

        # Verbose
        if verbose:
            print('bjobs (status):', ran.stdout)

        # Initialize job status
        job_status = 'EXIT'
        # Split bjobs row according to whitespaces
        row = re.split(r'\s+', ran.stdout)
        # Case row format match expected
        if len(row) >= 3:
            # Get job status
            job_status = row[2]

        # Return retrieved status
        return job_status

    @classmethod
    def kill(cls, job_id, verbose=False):
        # Retrieve command output
        ran = subprocess.run(
            check=True,
            capture_output=True,
            encoding='utf-8',
            args=['bkill', job_id]
        )
        # Verbose
        if verbose:
            print('Killed job {:s}'.format(str(job_id)))

    @classmethod
    def id_from_string(cls, in_string, verbose=False):
        """ Retrieve process id from verbose output

        Args
        in_string (str)     Input string, generated by submitting a job to lsf
        verbose (bool)      Whether to print out verbose log or not

        Return
        (str):              Job id, as string (avoids cutting leading zeroes)
        """
        # Initialize job id
        job_id = None
        # Split string according to newline
        in_list = in_string.split('\n')
        # Loop through each line
        for i in range(len(in_list)):
            # Case string matches job id
            match_id = re.search(r'^Job \<(\d+)\>', in_list[i])
            # Case there is a match
            if match_id:
                # Retrieve job id
                job_id = str(match_id.group(1))
                # Exit loop
                break

        # Verbose log
        if verbose:
            print('Retrieved job id {:s}'.format(str(job_id)))
            print('from {:s}'.format(in_string))

        # Return retrieved job id
        return job_id
