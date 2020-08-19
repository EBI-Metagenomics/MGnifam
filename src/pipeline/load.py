from src.scheduler import LSFScheduler, LocalScheduler
from src.mgnifam import Cluster
from src.hmm.hmmpress import HMMPress
from src.hmm.hmmemit import HMMEmit
from src.hmm.hmmer import Tblout, Domtblout
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
# from time import time
from glob import iglob
from time import time
import argparse
import os
import re


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
    def run(self, clusters_iter, mgnifam_version, verbose=False):
        """ Load clusters to MGnifam database

        First, iterate through each given cluster path and load it into MGnifam
        dataset. Then, chech MGnifam dataset for intersections among clusters.

        Note: either insertions and checks are done inside the same
        transaction. If some error happens, database is rollback to previous
        state.

        Args
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
                clusters_iter,
                verbose=verbose
            )

            # # Initialize sequence and domain scores dictionary
            # sequence_scores = dict()
            # domain_scores = dict()
            # # Loop through each cluster
            # for cluster_name in clusters.keys():
            #     # Get current cluster
            #     cluster = clusters[cluster_name]
            #     # Update sequence scores dictionary
            #     sequence_scores[cluster_name] = cluster.sequence_scores
            #     # Update domain scores dictionary
            #     domain_scores[cluster_name] = cluster.domain_scores
            #
            #
            #
            # # Get all HMM accession from MGnifam
            # mgnifam_acc = self.mgnifam_db.get_accessions()
            #
            # # Define MGnifam HMM path
            # mgnifam_hmm_path = os.path.join(clusters_path, 'MGnifam.hmm')
            # # Fill MGnifam HMM library
            # mgnifam_hmm = self.mgnifam_db.make_hmmlib(
            #     clusters_acc=mgnifam_acc,
            #     hmm_path=mgnifam_hmm_path
            # )
            #
            # # Define MGnifam consensus path
            # mgnifam_cns_path = os.path.join(clusters_path, 'MGnifam.fasta')
            # # Make MGnifam consensus
            # self.make_consensus(
            #     clusters_acc=mgnifam_acc,
            #     clusters_path=clusters_path,
            #     fasta_path=mgnifam_cns_path
            # )
            #
            # # Prepare HMM library for scan
            # self.prepare_scan(
            #     hmm_path=mgnifam_hmm.path,
            #     verbose=verbose
            # )
            #
            # # Define path to scan ouptut files
            # mgnifam_out_path = os.path.join(clusters_path, 'MGnifam.out')
            # mgnifam_seq_path = os.path.join(clusters_path, 'MGnifam.tblout')
            # mgnifam_dom_path = os.path.join(clusters_path, 'MGnifam.domtblout')
            # # Run mgnifam scan against consensus sequence
            # self.make_scan(
            #     # Set path to input HMM model/library
            #     hmm_path=mgnifam_hmm.path,
            #     # Set path to input FASTA file
            #     fasta_path=os.path.join(clusters_path)
            # )

        # If any exception got catched, rollback
        except Exception:
            # Discard changes
            self.mgnifam_db.rollback()
            # Go on raising exception
            raise

    # def get_clusters(self, clusters_path, verbose=False):
    #     """ Initialize clusters dictionary
    #
    #     Loop through each cluster directory in given clusters path.
    #     Moreover, it checks for foundamental files: HMM model file, domain
    #     hits against MGnifam and SEED alignment and ALIGN full alignment.
    #
    #     Args
    #     clusters_path (str)         Path where clusters DESC files are stored
    #     verbose (bool)              Wether to print verbose output
    #
    #     Return
    #     (dict)                      Dictionary mapping cluster name to cluster
    #                                 instance
    #     """
    #     # Initialize clusters dictionary (accession: cluster)
    #     clusters = dict()
    #     # Loop through each DESC file in clusters directories
    #     for cluster_path in iglob(os.path.join(clusters_path, 'MGYP*')):
    #         # Load cluster from DESC file
    #         cluster = Cluster.from_dir(cluster_path)
    #
    #         # Define HMM model path
    #         model_path = os.path.join(cluster_path, 'HMM.model')
    #
    #         # Define SEED alignment path
    #         seed_path = os.path.join(cluster_path, 'SEED.aln')
    #
    #         # Define ALIGN alignment path
    #         align_path = os.path.join(cluster_path, 'ALIGN.aln')
    #
    #         # Store cluster
    #         clusters[cluster.name] = cluster
    #
    #
    #     # Verbose output
    #     if verbose:
    #         print('Loaded {:d} clusters'.format(len(clusters)), end=' ')
    #         print('from {:s}'.fromat(clusters_path))
    #     # Return clusters dictionary
    #     return clusters

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
            print('clusters to MGnifam database...', end=' ')

        # Loop through each cluster path
        for cluster_path in clusters_iter:


            # Try loading cluster into database
            try:

                # Parse cluster from directory
                cluster = Cluster.from_dir(cluster_path)
                # Load cluster to MGnifam database
                cluster.to_mgnifam(self.mgnifam_db)

            # Case value error has been found
            except ValueError:
                # Go to next iteration
                pass

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Initialize cluster accession and id
        cluster_acc, cluster_id = '', ''
        # Loop through each cluster
        for cluster_name, cluster in clusters.items():

            # Case cluster accession is set: get next one
            if cluster_acc:
                # Remove prefix
                cluster_acc = int(re.sub(r'^MGYF', '', cluster_acc))
                # Increase accession number, reset prefix
                cluster_acc = 'MGYF{:05d}'.format(cluster_acc + 1)

            # Case cluster id is set: get next one
            if cluster_id:
                # Remove prefix, cast to integer
                cluster_id = int(re.sub(r'^MGDUF', '', cluster_id))
                # Increase id, reset prefix
                cluster_id = 'MGDUF{:04d}'.format(cluster_id + 1)

            # Update cluster accession number and id
            # NOTE: Avoid asking for new accession number and id again!
            cluster.acc = cluster_acc
            cluster.id = cluster_id

            # Define cluster SE
            cluster.seed = '{:s} ({:s})'.format(cluster_name, mgnifam_version)

            # Load cluster into database (updates cluster inplace)
            self.mgnifam_db.load_cluster(cluster)

            # TODO Load HMM model
            self.mgnifam_db.load_model(cluster.get_path('HMM'))

            # TODO Load SEED alignment
            self.mgnifam_db.load_seed(cluster.get_path('SEED'))

            # TODO Load ALIGN alignment
            self.mgnifam_db.load_align(cluster.get_path('ALIGN'))

            # TODO Load SEED alignment regions

            # TODO Load ALIGN alignment regions

            # Update current accession number and id
            cluster_acc = cluster.acc
            cluster_id = cluster.id

            # Verbose
            if verbose:
                print(' '.join([
                    'Cluster {:s}'.format(cluster.name),
                    'has been loaded into database'
                ]))

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            print(' '.join([
                'All {:d} clusters'.format(len(clusters)),
                'have been loaded',
                'in {:.0f} seconds'.format(time_end - time_beg)
            ]))

    def make_consensus(self, clusters_acc, fasta_path, verbose=False):
        """ Generate consensus FASTA file

        First, loop through each HMM in mgnifam, create HMM, run hmmemit
        Then, append fasta sequences retireved from hmmemit to consensus
        fasta file

        Args
        clusters_acc (set)          Set of cluster accession numbers
        fasta_path (str)            Path to output fasta file
        verbose (bool)              Whether to print verbose output
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Open consensus file in write mode
        consensus_file = open(fasta_path, 'w')
        # Loop through each cluster accession number in MGnifam
        for accession in clusters_acc:
            # Make new temporary HMM model file
            hmm_path = NamedTemporaryFile(delete=False, suffix='.hmm').name
            # Retrieve HMM for current cluster accession number
            self.mgnifam_db.make_hmm(
                accession=accession,
                path=hmm_path
            )

            # Define temporary output file
            out_path = NamedTemporaryFile(delete=False, suffix='.fasta').name
            # Run hmmemit for current cluster
            self.hmmemit(
                consensus=True,
                num_sequences=10,
                out_path=out_path,
                hmm_path=hmm_path
            )

            # Open temporary output file
            with open(out_path, 'r') as out_file:
                # Loop through each line in output file
                for line in out_file:
                    # Check if line is fasta header
                    is_header = re.search('^>', line)
                    # Case line is header
                    if is_header:
                        # Substitute line with current cluster accession
                        line = '>' + accession + '\n'
                    # Write line to consensus_file
                    consensus_file.write(line)

            # Remove temporary output file
            os.unlink(out_path)
            # Remove temporary HMM file
            os.unlink(hmm_path)

        # Update timers
        time_end = time()
        # Verbose
        if verbose:
            print(' '.join([
                'Done consensus fasta file',
                'for {:d} clusters'.format(len(clusters_acc)),
                'in {:.0f} seconds'.format(time_end - time_beg)
            ]))

        # Close consensus file
        consensus_file.close()

    def prepare_scan(self, hmm_path, verbose=False):
        """ Prepare files for hmmscan

        Args
        hmm_path (str)      Path where input HMM model/library is stored

        Return
        (str)               Path to generated .h3p file
        (str)               Path to generated .h3m file
        (str)               Path to generated .h3i file
        (str)               Path to generated .h3f file
        """
        # Try executing hmmpress
        try:
            # Run hmmpress script
            self.hmmpress(
                hmm_path=hmm_path,
                overwrite=True
            )

            # Get current directory
            hmm_dir = os.path.dirname(hmm_path)
            # Initialize generated file paths
            h3p_path = next(iglob(hmm_dir + '/*.h3p'))
            h3m_path = next(iglob(hmm_dir + '/*.h3m'))
            h3i_path = next(iglob(hmm_dir + '/*.h3i'))
            h3f_path = next(iglob(hmm_dir + '/*.h3f'))

            # Check if all the files are in place
            if not (h3p_path and h3m_path and h3i_path and h3f_path):
                # Raise file not found exception
                raise FileNotFoundError('\n'.join([
                    'Error: some files have not been generated:',
                    '- h3p: {:s}'.format(os.path.basename(h3p_path)),
                    '- h3m: {:s}'.format(os.path.basename(h3m_path)),
                    '- h3i: {:s}'.format(os.path.basename(h3i_path)),
                    '- h3f: {:s}'.format(os.path.basename(h3f_path))
                ]))

            # Return generated files path
            return h3p_path, h3m_path, h3i_path, h3f_path

        # Intercept process error
        except CalledProcessError as err:
            # Raise an explanatory exception
            raise Exception('\n'.join([
                'Error: hmmpress script returned {}:'.format(err.returncode),
                err.stderr.strip()
            ]))

    def make_scan(self, hmm_path, fasta_path, verbose=False):
        """ Scan MGnifam

        Scan all MGnifam clusters in root directory, using previously
        generated consensus sequences FASTA file and previously generated
        HMM model/library.

        Args
        hmm_path (str)          Path to input HMM model/library
        fasta_path (str)        Path to input FASTA file
        verbose (bool)          Whether to print verbose output

        Return
        (str)                   Path to output file
        (str)                   Path to tblout file
        (str)                   Path to domtblout file
        """
        # Try executing hmmscan
        try:

            # Get current HMM model directory
            hmm_dir = os.path.dirname(hmm_path)
            # Define output file paths
            out_path = hmm_dir + '/MGnifam.out'  # Output file
            tblout_path = hmm_dir + '/MGnifam.tblout'  # Sequence hits
            domtblout_path = hmm_dir + '/MGnifam.domtblout'  # Domain hits

            # Run hmmscan, make .domtblout and .tblout files
            self.hmmscan(
                # Set input HMM model path
                # NOTE: .h3[mifp] files must be in the same directory
                hmm_path=hmm_path,
                # Set input FASTA file path
                fasta_path=fasta_path,
                # Define output path
                out_path=out_path,
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

            # Return result paths
            return out_path, tblout_path, domtblout_path

        # Intercept subprocess error
        except CalledProcessError as err:
            # Raise new exception
            raise Exception('\n'.join([
                'Error: hmmscan returned {}'.format(err.returncode),
                err.stderr.strip()
            ]))

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
    # Define path to input directory
    parser.add_argument(
        '-i', '--in_path', nargs='+', type=str, required=True,
        help='Path to input cluster directories'
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
    # # Define E-value significance threshold for both UniProt and MGnifam
    # parser.add_argument(
    #     '-e', '--e_value', type=float, default=0.01,
    #     help='E-value threhsold for both UniProt and MGnifam comparisons'
    # )
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

    #
