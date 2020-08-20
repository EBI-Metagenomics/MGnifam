# from src.scheduler import LSFScheduler, LocalScheduler
from src.hmm.hmmer import Tblout, Domtblout
from src.hmm.hmmpress import HMMPress
from src.hmm.hmmemit import HMMEmit
from src.sequences import fasta_iter
from src.mgnifam import Cluster

from subprocess import CalledProcessError
from tempfile import mkstemp
from glob import glob, iglob
from time import time
import argparse
import shutil
import os


class Load(object):
    """ Load MGnifam relese into database

    Starting from a build directory, loads each new MGnifam cluster and
    searches for intersections with either previous MGnifam and Pfam release.
    """

    # Constructor
    def __init__(
        self, mgnifam_db, pfam_db,
        hmmpress_cmd=['hmmpress'],
        hmmemit_cmd=['hmmemit'],
        env=os.environ
    ):
        # Store databases
        self.mgnifam_db = mgnifam_db
        self.pfam_db = pfam_db
        # Define scripts instances
        self.hmmpress = HMMPress(cmd=hmmpress_cmd, env=env)
        self.hmmemit = HMMEmit(cmd=hmmemit_cmd, env=env)
        # Store environmental variables
        self.env = env

    # Run pipeline
    def run(self, release_path, clusters_iter, mgnifam_version, verbose=False):
        """ Load clusters to MGnifam database

        First, iterate through each given cluster path and load it into MGnifam
        dataset. Then, chech MGnifam dataset for intersections among clusters.

        Note: either insertions and checks are done inside the same
        transaction. If some error happens, database is rollback to previous
        state.

        Args
        release_path (str)          Path to release directory
        clusters_iter (list)        List of cluster directories paths
        mgnifam_version (str)       String for mgnifam version
        verbose (bool)              Whether to print out verbose log
        """

        # Try to make new release
        try:

            # Connect to database
            self.mgnifam_db.connect()

            # Load each cluster in database
            self.load_clusters(
                clusters_iter=clusters_iter,
                verbose=verbose
            )

            # Initialize release directory
            self.init_dir(
                release_path=release_path
            )

            # Retrieve MGnifam clusters
            mgnifam_dict = self.retrieve_clusters(
                release_path=release_path,
                verbose=verbose
            )
            # TODO Retrieve MGnifam clusters domain scores
            # TODO Retrieve MGnifam clusters sequence scores

            # Make MGnifam consensus
            self.make_consensus(
                clusters_path=release_path,
                verbose=verbose
            )

            # Make MGnifam HMM model
            self.prepare_scan(
                clusters_path=release_path,
                verbose=verbose
            )

            # Make MGnifam scan
            tblout_path, domtblout_path = self.make_scan(
                clusters_path=release_path,
                verbose=verbose
            )

            # TODO Parse MGnifam scan
            # TODO Kill clusters intercepting

            # Raise exception
            raise Exception()

        # If any exception got catched, rollback
        except Exception:
            # Discard changes
            self.mgnifam_db.rollback()
            # Go on raising exception
            raise

        # Otherwise, commit changes
        self.mgnifam_db.conn.commit()

    def load_clusters(self, clusters_iter, mgnifam_version, verbose=False):
        """ Load clusters into MGnifam database

        Args
        clusters (dict)     Dict of clusters retrieved from directory
        verbose (bool)      Whether to print verbose output
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Loading {:d}'.format(len(clusters_iter)), end=' ')
            print('clusters to MGnifam database...')

        # Loop through each cluster path
        for cluster_path in clusters_iter:
            # Try loading cluster into database
            try:
                # Parse cluster from directory
                cluster = Cluster.from_dir(cluster_path)
                # Load cluster to MGnifam database
                cluster.to_mgnifam(self.mgnifam_db)

            # Case value error has been found
            except ValueError as err:
                # Verbose
                if verbose:
                    # Print discarded cluster
                    print('Cluster {:s} discarded:'.format(cluster.name))
                    print(err.message.strip())
                # Go to next iteration
                pass

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('...done in {:.0f} seconds'.format(time_end - time_beg))

    def init_release(clusters_path):
        """ Initialize release directory structure

        Args
        clusters_path (str)         Path to release directory
        """
        # Check if release directory is already in place
        if not os.path.isdir(clusters_path):
            # Try making release directory
            os.mkdir(clusters_path)

    def retrieve_clusters(self, clusters_path, verbose=False):
        """ Retrieve all MGnifam clusters

        Loop through each accession in MGnifam dataset, retrieve it and make
        cluster directory. Fill directory with HMM.model file.

        Args
        clusters_path (str)         Path to release directory
        verbose (bool)              Whether to print out verbose log

        Return
        (dict)                      Dictionary mapping cluster names to
                                    clusters objects
        """
        # Set database reference
        db = self.mgnifam_db

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Retrieving MGnifam clusters...', end=' ')

        # Initialize clusters dictionary
        clusters_dict = dict()
        # Retrieve all sequence accessions from MGnifam database
        clusters_acc = db.get_accessions()
        # Loop through all retrieved accessions
        for acc in clusters_acc:
            # Define cluster path
            cluster_path = os.path.join(clusters_path, acc)

            # Try making cluster directory
            try:
                # Case cluster  path does not exist
                if not os.path.isdir(cluster_path):
                    # Make cluster directory
                    os.mkdir(cluster_path)

                # Retrieve cluster
                cluster = db.get_cluster(acccession=acc)
                # Set cluster path
                cluster.path = cluster_path
                # Make DESC file
                cluster.to_desc()

                # Define HMM model path
                model_path = os.path.join(cluster_path, 'HMM.model')
                # Make HMM model
                db.make_hmm_model(accession=acc, path=model_path)

                # Append cluster to clusters dictionary
                clusters_dict.setdefault(acc, cluster)

            # Catch exceptions
            except Exception:
                # Discard cluster directory
                os.remove(cluster_path)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:0f} seconds'.format(time_end -  time_beg))

        # Return clusters dictionary
        return clusters_dict

    def make_consensus(self, clusters_path, verbose=False):
        """ Generate consensus FASTA file

        First, loop through each HMM in mgnifam and run hmmemit
        Then, append fasta sequences retireved from hmmemit to consensus
        fasta file

        Args
        clusters_path (set)         Path to folder holding clusters directories
        verbose (bool)              Whether to print verbose output
        """
        # Define HMM models iterator
        models_iter = glob(os.path.join(clusters_path, 'MGYF*', 'HMM.model'))

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making consensus for', edn=' ')
            print('{:d} clusters'.format(len(models_iter)), end=' ')
            print('at {:s}...'.format(clusters_path), end=' ')

        # Define all clusters FASTA consensus file
        cons_path = os.path.join(clusters_path, 'CONS.fa')
        # Open all-clusters consensus FASTA file
        with open(cons_path, 'w') as cons_file:
            # Loop through each HMM model
            for model_path in models_iter:
                # Get cluster path
                cluster_path = os.path.dirname(model_path)
                # Define consensus fasta file path
                fasta_path = os.path.join(cluster_path, 'CONS.fa')

                # Just call HMM emit
                self.hmmemit(
                    model_path=model_path,
                    out_path=fasta_path
                )

                # Define consensus fasta sequences
                fasta_sequences = dict()
                # Open output fasta file
                with open(fasta_path, 'r') as fasta_file:
                    # Iterate through each fasta entry in given file
                    for fasta_entry in fasta_iter(fasta_file):
                        # Split entry into header and residues
                        head, resid = tuple(fasta_entry.split('\n'))
                        # Store fasta sequence
                        fasta_sequences.setdefault(head, resid)

                # Get cluster accession (cluster folder is named after accession)
                cluster_acc = os.path.basename(cluster_path)
                # Open fasta file in write mode
                with open(fasta_path, 'w') as fasta_file:
                    # Iterate through each stored fasta tuple(header, residues)
                    for head, resid in fasta_sequences.items():
                        # Write out fasta line(s)
                        fasta_file.write('>' + cluster_acc + '\n' + resid + '\n')

                # Open fasta file in read mode
                with open(fasta_path, 'r') as fasta_file:
                    # Copy model to library
                    shutil.copyfileobj(fasta_file, cons_file)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    def prepare_scan(self, clusters_path, verbose=False):
        """ Prepare files for hmmscan

        First, make single HMM library by joining together HMM models stored
        in each cluster folder. Then, generate `.h3[pmif]` HMM dataset files.

        Args
        clusters_path (set)         Path to folder holding clusters directories
        verbose (bool)              Whether to print verbose output

        Raise
        (FileNotFoundError)         In case some output file has not been found
        (Exception)                 In case generic exception happened, such as
                                    errors during hmmpress execution
        """
        # Define clusters HMM model iterator
        models_iter = glob(os.path.join(clusters_path, 'MGYF*', 'HMM.model'))

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making HMM library and database', end=' ')
            print('out of {:d} HMM models'.format(len(models_iter)), end=' ')
            print('at {:s}...'.format(clusters_path), end=' ')

        # Define HMM library path
        library_path = os.path.join(clusters_path, 'HMM.library')
        # Open output file
        with open(library_path, 'w') as library_file:
            # Loop through each HMM model file
            for model_path in models_iter:
                # Open input file
                with open(model_path, 'r') as model_file:
                    # Copy model to library
                    shutil.copyfileobj(model_file, library_file)

        # Try running hmmpress
        try:
            # Run hmmpress script against library
            self.hmmpress(model_path=library_path, overwrite=True)
            # Initialize generated file paths
            h3p_path = next(iglob(clusters_path + '/*.h3p'))
            h3m_path = next(iglob(clusters_path + '/*.h3m'))
            h3i_path = next(iglob(clusters_path + '/*.h3i'))
            h3f_path = next(iglob(clusters_path + '/*.h3f'))

            # Check if all the files are in place
            if not (h3p_path and h3m_path and h3i_path and h3f_path):
                # Raise file not found exception
                raise FileNotFoundError('\n'.join([
                    'some files have not been generated:',
                    '- h3p: {:s}'.format(os.path.basename(h3p_path)),
                    '- h3m: {:s}'.format(os.path.basename(h3m_path)),
                    '- h3i: {:s}'.format(os.path.basename(h3i_path)),
                    '- h3f: {:s}'.format(os.path.basename(h3f_path))
                ]))

        # Intercept process error
        except CalledProcessError as err:
            # Raise exception
            raise Exception('\n'.join([
                'hmmpress script returned {}:'.format(err.returncode),
                err.stderr.strip()
            ]))

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    def make_scan(self, clusters_path, verbose=False):
        """ Scan MGnifam

        Scan all MGnifam clusters in root directory, using previously
        generated consensus sequences FASTA file and previously generated
        HMM model/library.

        Args
        clusters_path (str)         Path to clusters directory
        verbose (bool)              Whether to print verbose output

        Return
        (str)                       Path to output file
        (str)                       Path to tblout file
        (str)                       Path to domtblout file

        Raise
        (FileNotFoundError)         In case some file has not been found
        (Exception)                 In case generic exception happened, such as
                                    errors during hmmpress execution
        """
        # Define HMM library path
        library_path = os.path.join(clusters_path, 'HMM.library')
        # Initialize FASTA consensus file path
        cons_path = os.path.join(clusters_path, 'CONS.fa')

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Scanning clusters at', end=' ')
            print('{:s}...'.format(clusters_path), end=' ')

        # Case cannot find HMM library file
        if not os.path.isfile(library_path):
            # Raise exception
            raise FileNotFoundError(
                'could not find HMM library at {:s}'.format(library_path)
            )

        # Case cannot find HMM library file
        if not os.path.isfile(cons_path):
            # Raise exception
            raise FileNotFoundError(
                'could not find consensus FASTA file at {:s}'.format(cons_path)
            )

        # Try executing hmmscan
        try:
            # Define tblout path
            tblout_path = mkstemp(dir=clusters_path, suffix='.tblout')
            # Define domtblout path
            domtblout_path = mkstemp(dir=clusters_path, suffix='.domtblout')
            # Run hmmscan, make .domtblout and .tblout files
            self.hmmscan(
                # Set input HMM model path
                # NOTE: .h3[mifp] files must be in the same directory
                hmm_path=library_path,
                # Set input FASTA file path
                fasta_path=cons_path,
                # # Define output path
                # out_path=out_path,
                # Define sequence hits output path
                tblout_path=tblout_path,
                # Define domain hits output path
                domtblout_path=domtblout_path,
                # Do not compute alignments
                no_alignments=True,
                # Set an high reporting threshold for sequences
                sequenece_eval=(1000, ),
                # Set an high reporting threshold for domains
                domain_eval=(1000, )
            )

        # Intercept hmmscan execution error
        except CalledProcessError as err:
            # Raise exception
            raise Exception('hmmscan returned {}:\n{}'.format(
                                err.returncode,
                                err.stderr.strip()
                            ))

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Return either tblout and domtblout paths
        return tblout_path, domtblout_path

    def parse_scan(self, sequence_scores, domain_scores, tblout_path, domtblout_path, verbose=False):
        """ Parse hmmscan output files

        Args
        sequence_scores (dict)  Dictionary mappung cluster accession numbers
                                to sequences (TC, NC, GA) scores
        domain_scores (dict)    Dictionary mapping cluster accession numbers
                                to domains (TC, NC, GA) scores
        tblout_path (str)       Path to sequence hits file
        domtblout_path (str)    Path to domain hits file
        verbose (bool)          Whether to print verbose output

        Return
        (dict)                  Dictionary mapping cluster accession numbers
                                (keys) to set of other intersecting cluster
                                accession numbers (values)
        """
        # Parse tblout
        tblout = Tblout(tblout_path)
        # Retireve hits
        sequence_hits = tblout.get_hits(scores=sequence_scores)

        # Parse domtblout
        domtblout = Domtblout(domtblout_path)
        # Retrieve hits
        domain_hits = domtblout.get_hits(scores=domain_scores)

        # # Keep only hits in both sequences and domains
        # query_names = set(sequence_hits.keys()) & set(domain_hits.keys())
        # # Loop through each query name
        # for query_name in query_names:
        #     # Retrieve each domain hit
        #     for
        raise NotImplementedError


# Command line execution
if __name__ == '__main__':

    # Argparse
    parser = argparse.ArgumentParser(description='Build MGnifam clusters')
    # Define path to input clusters files
    parser.add_argument(
        '-i', '--in_path', nargs='+', type=str, required=True,
        help='Path to input cluster directories'
    )
    # Define path to release directory
    parser.add_argument(
        '-o', '--out_path', type=str, required=True,
        help='Path to release directory'
    )
    # # Define output path
    # parser.add_argument(
    #     '-o', '--out_path', type=str, required=True,
    #     help='Path to output directory'
    # )
    # # Define author name
    # parser.add_argument(
    #     '-a', '--author_name', type=str, default='Anonymous',
    #     help='Name of the user running MGnifam build pipeline'
    # )
    # # Wether to shuffle input or not
    # parser.add_argument(
    #     '--shuffle', type=int, default=0,
    #     help='Whether to shuffle input cluster names or not'
    # )
    # # Define batch size
    # parser.add_argument(
    #     '--batch_size', type=int, default=100,
    #     help='Number of clusters to make at each iteration'
    # )
    # # Define maximum number of clusters to make
    # parser.add_argument(
    #     '--max_clusters', type=int, default=0,
    #     help='Maximum number of clusters to make'
    # )
    # # Define path to LinClust clusters file(s)
    # parser.add_argument(
    #     '--linclust_path', nargs='+', type=str, default=glob(LINCLUST_PATH),
    #     help='Path to LinClust clusters file(s)'
    # )
    # # Define path to MGnifam file(s)
    # parser.add_argument(
    #     '--mgnifam_path', nargs='+', type=str, default=glob(MGNIFAM_PATH),
    #     help='Path to MGnifam file(s)'
    # )
    # # Define MGnifam width (maximum sequence length)
    # parser.add_argument(
    #     '--mgnifam_width', type=int, required=False,
    #     help='Maximum sequence length in MGnifam dataset'
    # )
    # # Define MGnifam height (total number of sequences)
    # parser.add_argument(
    #     '--mgnifam_height', type=int, required=False,
    #     help='Total number of sequences in MGnifam dataset'
    # )
    # # Define path to UniProt file(s)
    # parser.add_argument(
    #     '--uniprot_path', nargs='+', type=str, default=glob(UNIPROT_PATH),
    #     help='Path to UniProt file(s)'
    # )
    # # Define UniProt width (maximum sequence length)
    # parser.add_argument(
    #     '--uniprot_width', type=int, required=False,
    #     help='Maximum sequence length in UniProt dataset'
    # )
    # # Define UniProt height (total number of sequences)
    # parser.add_argument(
    #     '--uniprot_height', type=int, required=False,
    #     help='Total number of sequences in UniProt dataset'
    # )
    # # Define MobiDB executable
    # parser.add_argument(
    #     '--mobidb_cmd', nargs='+', type=str, default=['mobidb_lite.py'],
    #     help='MobiDB Lite disorder predictor executable'
    # )
    # # Define Muscle executable
    # parser.add_argument(
    #     '--muscle_cmd', nargs='+', type=str, default=['muscle'],
    #     help='Muscle multiple sequence aligner executable'
    # )
    # # Define hmmsearch executable
    # parser.add_argument(
    #     '--hmmsearch_cmd', nargs='+', type=str, default=['hmmsearch'],
    #     help='HMMER3 search executable'
    # )
    # # Define hmmbuild executable
    # parser.add_argument(
    #     '--hmmbuild_cmd', nargs='+', type=str, default=['hmmbuild'],
    #     help='HMMER3 build executable'
    # )
    # # Define hmmalign executable
    # parser.add_argument(
    #     '--hmmalign_cmd', nargs='+', type=str, default=['hmmalign'],
    #     help='HMMER3 align executable'
    # )
    # Define hmmpress executable
    parser.add_argument(
        '--hmmpress_cmd', nargs='+', type=str, default=['hmmpress'],
        help='HMMER3 press executable'
    )
    # Define hmmemit executable
    parser.add_argument(
        '--hmmemit_cmd', nargs='+', type=str, default=['hmmemit'],
        help='HMMER3 emit executable'
    )
    # Whether to print verbose output
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Print verbose output'
    )
    # Define E-value significance threshold for both MGnifam and Pfam
    parser.add_argument(
        '-e', '--e_value', type=float, default=0.01,
        help='E-value threhsold for both UniProt and MGnifam comparisons'
    )
    # Define environmental variables .json file
    parser.add_argument(
        '--env_path', type=str, required=False,
        help='Path to .json file holding environmental variables'
    )
    # Define MGnifam database options
    group = parser.add_argument_group('MGnifam database options')
    # Add database user
    group.add_argument(
        '--mgnifam_user', type=str, required=True,
        help='MGnifam database username'
    )
    # Add database password
    group.add_argument(
        '--mgnifam_password', type=str, required=False,
        help='MGnifam database password'
    )
    # Add database host
    group.add_argument(
        '--mgnifam_host', type=str, default='localhost',
        help='MGnifam database host'
    )
    # Add database host
    group.add_argument(
        '--mgnifam_port', type=str, default=3309,
        help='MGnifam database host'
    )
    # Add database user
    group.add_argument(
        '--pfam_user', type=str, required=True,
        help='Pfam database username'
    )
    # Add database password
    group.add_argument(
        '--pfam_password', type=str, required=False,
        help='Pfam database password'
    )
    # Add database host
    group.add_argument(
        '--pfam_host', type=str, default='localhost',
        help='Pfam database host'
    )
    # Add database host
    group.add_argument(
        '--pfam_port', type=str, default=3309,
        help='Pfam database host'
    )
    # # Define scheduler options
    # group = parser.add_argument_group('Scheduler options')
    # # Define schduler type
    # group.add_argument(
    #     '-s', '--scheduler_type', type=str, default='LSF',
    #     help='Type of scheduler to use to distribute parallel processes'
    # )
    # # Define minimum number of jobs
    # group.add_argument(
    #     '-j', '--min_jobs', type=int, default=0,
    #     help='Minimum number of parallel processes to keep alive'
    # )
    # # Define maximum number of jobs
    # group.add_argument(
    #     '-J', '--max_jobs', type=int, default=100,
    #     help='Maximum number of parallel processes to keep alive'
    # )
    # # Define minimum number of cores
    # group.add_argument(
    #     '-c', '--min_cores', type=int, default=1,
    #     help='Minimum number of cores to use per process'
    # )
    # # Define maximum number of cores
    # group.add_argument(
    #     '-C', '--max_cores', type=int, default=1,
    #     help='Maximum number of cores to use per process'
    # )
    # # Define minimum memory allocable per job
    # group.add_argument(
    #     '-m', '--min_memory', type=str, default='2 GB',
    #     help='Minimum memory allocable per process'
    # )
    # # Define maximum memory allocable per job
    # group.add_argument(
    #     '-M', '--max_memory', type=str, default='4 GB',
    #     help='Maximum memory allocable per process'
    # )
    # # Define walltime
    # group.add_argument(
    #     '-W', '--walltime', type=str, default='02:00',
    #     help='How long can a process be kept alive'
    # )
    # Retrieve arguments
    args = parser.parse_args()

    # # Initialize scheduler
    # scheduler = None
    # # Case scheduler type is LSF
    # if args.scheduler_type == 'LSF':
    #     # Define LSF scheduler
    #     scheduler = LSFScheduler(
    #         # Define cores boundaries
    #         min_cores=args.min_cores,
    #         max_cores=args.max_cores,
    #         # Define memory boundaries
    #         min_memory=args.min_memory,
    #         max_memory=args.max_memory,
    #         # Debug
    #         silence_logs='debug',
    #         # Define walltime
    #         walltime=args.walltime,
    #         # Define processes per job
    #         processes=1
    #     )
    # # Case scheduler type is Local
    # if args.scheduler_type == 'Local':
    #     # Define Local scheduler
    #     scheduler = LocalScheduler(
    #         # Define cores boundaries
    #         min_cores=args.min_cores,
    #         max_cores=args.max_cores,
    #         # Define memory boundaries
    #         min_memory=args.min_memory,
    #         max_memory=args.max_memory,
    #         # Define processes per job
    #         threads_per_worker=1
    #     )
    # # Case no scheduler has been set
    # if scheduler is None:
    #     # Raise new error
    #     raise ValueError(' '.join([
    #         'scheduler type can be one among `Local` or `LSF`',
    #         '%s has been chosen instead' % args.scheduler_type
    #     ]))
