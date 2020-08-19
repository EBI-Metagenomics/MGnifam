# Dependencies
from src.scheduler import LSFScheduler, LocalScheduler
from src.pipeline.pipeline import Pipeline
from src.dataset import LinClust, Fasta
from time import time
import argparse
import os


class Dataset(Pipeline):
    """Initialize dataset

    Starting from full, big size datasets, splits them in chunks.
    Actually considered dataset are MGnifam, UniProt and LinClust clusters.
    """

    # Constructor
    def __init__(self, scheduler):
        # Save scheduler instance
        self.scheduler = scheduler

    # Execute pipeline
    def run(
        # Number of jobs, verbose log
        self, min_jobs, max_jobs, verbose=False,
        # LinClust dataset parameters
        linclust_in_path='', linclust_chunk_size=1e07, linclust_out_path='{:04d}',
        # UniProt dataset parameters
        uniprot_in_path='', uniprot_chunk_size=1e06, uniprot_out_path='{:04d}',
        # MGnify dataset parameters
        mgnifam_in_path='', mgnifam_chunk_size=1e06, mgnifam_out_path='{:04d}'
    ):

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making dataset chunks...')

        # Adapt scheduler
        self.scheduler.adapt(min_jobs, max_jobs)
        # Open scheduler connection to client
        with self.scheduler.client as client:

            # Initialize futures dict(dataset_name: future)
            futures = dict()

            # Case linclust path is defined
            if linclust_in_path:
                # Define single chunk LinClust dataset
                linclust = LinClust(path=linclust_in_path)
                # Split in chunks
                futures['linclust'] = client.submit(
                    # Chunking method
                    linclust.to_chunks,
                    # Chunks output path,
                    chunk_path=linclust_out_path,
                    # Size (number of lines) of every cluster
                    chunk_size=linclust_chunk_size
                )
                # Verbose
                if verbose:
                    print('  chunking {:s}'.format(linclust_in_path), end=' ')
                    print('dataset in {:d} chunks'.format(linclust_chunk_size))

            # Case uniprot path is defined
            if uniprot_in_path:
                # Define single chunk UniProt dataset
                uniprot = Fasta(path=uniprot_in_path)
                # Split in chunks
                futures['uniprot'] = client.submit(
                    # Chunking method
                    uniprot.to_chunks,
                    # Chunks output path,
                    chunk_path=uniprot_out_path,
                    # Size (number of fasta entries) of every cluster
                    chunk_size=uniprot_chunk_size
                )
                # Verbose
                if verbose:
                    print('  chunking {:s}'.format(uniprot_in_path), end=' ')
                    print('dataset in {:d} chunks'.format(uniprot_chunk_size))

            # Case mgnify path is defined
            if mgnifam_in_path:
                # Define single chunk MGnifam dataset
                mgnifam = Fasta(path=mgnifam_in_path)
                # Split in chunks
                futures['mgnifam'] = client.submit(
                    # Chunking method
                    mgnifam.to_chunks,
                    # Chunks output path,
                    chunk_path=mgnifam_out_path,
                    # Size (number of fasta entries) of every cluster
                    chunk_size=mgnifam_chunk_size
                )
                # Verbose
                if verbose:
                    print('  chunking {:s}'.format(mgnifam_in_path), end=' ')
                    print('dataset in {:d} chunks'.format(mgnifam_chunk_size))

            # Wait for every future to finish
            client.gather(futures)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))


# Command line execution
if __name__ == '__main__':

    # Define project root path
    ROOT_PATH = os.path.dirname(__file__) + '/../..'
    # Define path to data folder
    DATA_PATH = ROOT_PATH + '/data'
    # Define path to LinClust file(s)
    LINCLUST_PATH = DATA_PATH + '/clusters/chunk*.tsv.gz'
    # Define path to MGnifam file(s)
    MGNIFAM_PATH = DATA_PATH + '/mgnify/chunk*.fa.gz'
    # Define path to UniProt file(s)
    UNIPROT_PATH = DATA_PATH + '/uniprot/chunk*.fa.gz'

    # Argparse
    parser = argparse.ArgumentParser(description='Prepare dataset')
    # Define path to LinClust clusters dataset
    parser.add_argument(
        '--linclust_path', type=str, required=False,
        help='Path to LinClust clusters file(s)'
    )
    # Define LinClust chunk size
    parser.add_argument(
        '--linclust_chunk_size', type=int, default=int(1e07),
        help='LinClust chunk size'
    )
    # Define path to MGnifam dataset
    parser.add_argument(
        '--mgnifam_path', type=str, required=False,
        help='Path to MGnifam file'
    )
    # Define MGnifam chunk size
    parser.add_argument(
        '--mgnifam_chunk_size', type=int, default=int(1e06),
        help='MGnifam chunk size'
    )
    # Define path to UniProt dataset
    parser.add_argument(
        '--uniprot_path', type=str, required=False,
        help='Path to UniProt file'
    )
    # Define UniProt chunk size
    parser.add_argument(
        '--uniprot_chunk_size', type=int, default=int(1e06),
        help='UniProt chunk size'
    )
    # Define output directory path
    parser.add_argument(
        '-o', '--out_path', type=str, default=DATA_PATH,
        help='Path to output directory'
    )
    # Wether to gunzip chunks or not
    parser.add_argument(
        '-gz', '--gzip', type=int, default=0,
        help='Whether to compress chunks or not'
    )
    # Whether to print verbose output
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Print verbose output'
    )
    # Define scheduler options
    group = parser.add_argument_group('Scheduler options')
    # Define schduler type
    group.add_argument(
        '-s', '--scheduler_type', type=str, default='LSF',
        help='Type of scheduler to use to distribute parallel processes'
    )
    # Define minimum number of jobs
    group.add_argument(
        '-j', '--min_jobs', type=int, default=0,
        help='Minimum number of parallel processes to keep alive'
    )
    # Define maximum number of jobs
    group.add_argument(
        '-J', '--max_jobs', type=int, default=100,
        help='Maximum number of parallel processes to keep alive'
    )
    # Define minimum number of cores
    group.add_argument(
        '-c', '--min_cores', type=int, default=1,
        help='Minimum number of cores to use per process'
    )
    # Define maximum number of cores
    group.add_argument(
        '-C', '--max_cores', type=int, default=1,
        help='Maximum number of cores to use per process'
    )
    # Define minimum memory allocable per job
    group.add_argument(
        '-m', '--min_memory', type=str, default='1 GB',
        help='Minimum memory allocable per process'
    )
    # Define maximum memory allocable per job
    group.add_argument(
        '-M', '--max_memory', type=str, default='1 GB',
        help='Maximum memory allocable per process'
    )
    # Define walltime
    group.add_argument(
        '-W', '--walltime', type=str, default='02:00',
        help='How long can a process be kept alive'
    )
    # Retrieve arguments
    args = parser.parse_args()

    # Initialize scheduler
    scheduler = None
    # Case scheduler type is LSF
    if args.scheduler_type == 'LSF':
        # Define LSF scheduler
        scheduler = LSFScheduler(
            # Define cores boundaries
            min_cores=args.min_cores,
            max_cores=args.max_cores,
            # Define memory boundaries
            min_memory=args.min_memory,
            max_memory=args.max_memory,
            # Debug
            silence_logs='debug',
            # Define walltime
            walltime=args.walltime,
            # Define processes per job
            processes=1
        )
    # Case scheduler type is Local
    if args.scheduler_type == 'Local':
        # Define Local scheduler
        scheduler = LocalScheduler(
            # Define cores boundaries
            min_cores=args.min_cores,
            max_cores=args.max_cores,
            # Define memory boundaries
            min_memory=args.min_memory,
            max_memory=args.max_memory,
            # Define processes per job
            threads_per_worker=1
        )
    # Case no scheduler has been set
    if scheduler is None:
        # Raise new error
        raise ValueError(' '.join([
            'scheduler type can be one among `Local` or `LSF`',
            '%s has been chosen instead' % args.scheduler_type
        ]))

    # Initialize pipeline
    pipeline = Dataset(scheduler=scheduler)

    # Initialize parameters dictionary
    params = dict()
    # Retrieve output directory path
    out_path = args.out_path

    # Case LinClust path is set
    if args.linclust_path:
        # Define LinClust output
        linclust_chunk = os.path.join(out_path, 'linclust', 'chunk{:06}.tsv')
        # Eventually add '.gz' extension to compress chunk
        linclust_chunk += '.gz' if args.gzip else ''
        # Make chunks directory
        os.makedirs(os.path.dirname(linclust_chunk), exist_ok=True)
        # Add LinClust parameters
        params = {**params, **{
            'linclust_in_path': args.linclust_path,
            'linclust_out_path': linclust_chunk,
            'linclust_chunk_size': args.linclust_chunk_size
        }}

    # Case UniProt path is set
    if args.uniprot_path:
        # Define UniProt chunk path
        uniprot_chunk = os.path.join(out_path, 'uniprot', 'chunk{:06}.fa')
        # Eventually add '.gz' extension to compress chunk
        uniprot_chunk += '.gz' if args.gzip else ''
        # Make chunks directory
        os.makedirs(os.path.dirname(uniprot_chunk), exist_ok=True)
        # Add UniProt parameters
        params = {**params, **{
            'uniprot_in_path': args.uniprot_path,
            'uniprot_out_path': uniprot_chunk,
            'uniprot_chunk_size': args.uniprot_chunk_size
        }}

    # Case MGnifam path is set
    if args.mgnifam_path:
        # Define MGnifam chunk path
        mgnifam_chunk = os.path.join(out_path, 'mgnifam', 'chunk{:06}.fa')
        # Eventually add '.gz' extension to compress chunk
        mgnifam_chunk += '.gz' if args.gzip else ''
        # Make chunks directory
        os.makedirs(os.path.dirname(mgnifam_chunk), exist_ok=True)
        # Add MGnifam parameters
        params = {**params, **{
            'mgnifam_in_path': args.mgnifam_path,
            'mgnifam_out_path': mgnifam_chunk,
            'mgnifam_chunk_size': args.mgnifam_chunk_size
        }}

    # Run the pipeline
    pipeline(
        **params,
        min_jobs=args.min_jobs,
        max_jobs=args.max_jobs,
        verbose=bool(args.verbose)
    )
