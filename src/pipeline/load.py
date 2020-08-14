from src.mgnifam import Cluster
from src.hmm.hmmpress import HMMPress
from src.hmm.hmmemit import HMMEmit
from src.hmm.hmmer import Tblout, Domtblout
from subprocess import CalledProcessError
from tempfile import NamedTemporaryFile
# from time import time
from glob import iglob
from time import time
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
        env=os.environ.copy()
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
    def run(self, clusters_path, mgnifam_version, verbose=False):

        # Try to make new release
        try:

            # Connect to database
            self.mgnifam_db.connect()

            # Initialize transaction
            self.mgnifam_db.autocommit = False

            # Retrieve clusters dictionary (accession: cluster)
            clusters = self.get_clusters(clusters_path, verbose=verbose)

            # Load each cluster in database
            self.load_clusters(clusters, verbose=verbose)

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

    def get_clusters(self, clusters_path, verbose=False):
        """ Initialize clusters dictionary

        Args
        clusters_path (str)         Path where clusters DESC files are stored
        verbose (bool)              Wether to print verbose output

        Return
        (dict)                      Dictionary mapping cluster name to cluster
                                    instance
        """
        # Initialize clusters dictionary (accession: cluster)
        clusters = dict()
        # Loop through each DESC file in clusters directories
        for cluster_path in iglob(os.path.join(clusters_path, 'MGYP*')):
            # Load cluster from DESC file
            cluster = Cluster.from_dir(cluster_path)
            # Store cluster
            clusters[cluster.name] = cluster
        # Verbose output
        if verbose:
            print('Loaded {:d} clusters'.format(len(clusters)), end=' ')
            print('from {:s}'.fromat(clusters_path))
        # Return clusters dictionary
        return clusters

    def load_clusters(self, clusters, mgnifam_version, verbose=False):
        """ Load clusters into MGnifam database

        Args
        clusters (dict)     Dict of clusters retrieved from directory
        verbose (bool)      Whether to print verbose output
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
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
