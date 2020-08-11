from src.mgnifam import Cluster
# from src.hmm.hmm import HMM
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
        # Store commands
        self.hmmpress_cmd = hmmpress_cmd,
        self.hmmemit_cmd = hmmemit_cmd,
        # Store environmental variables
        self.env = env

    # Run pipeline
    def run(self, clusters_path, verbose=False):

        # Retrieve clusters dictionary (accession: cluster)
        clusters = self.get_clusters(clusters_path, verbose=verbose)

        # Connect to database
        self.mgnifam_db.connect()

        # TODO Initialize transaction (locks database for a while)

        # Get all HMM accession from MGnifam
        mgnifam_acc = self.mgnifam_db.get_accessions()

        # Fill MGnifam HMM library
        mgnifam_hmm = self.mgnifam_db.make_hmmlib(
            clusters_acc=mgnifam_acc,
            hmm_path=os.path.join(clusters_path, 'MGNIFAM.hmm')
        )

        # Prepare HMM library for scan
        self.prepare_scan(
            hmm_path=mgnifam_hmm.path,
            verbose=verbose
        )

        # Make MGnifam consensus
        self.make_consensus(
            clusters_acc=mgnifam_acc,
            clusters_path=clusters_path
        )

        # TODO Run mgnifam scan

        raise NotImplementedError

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
        for desc_path in iglob(os.path.join(clusters_path, 'MGYP*', 'DESC')):
            # Load cluster from DESC file
            cluster = Cluster.from_desc(desc_path)
            # Store cluster
            clusters[cluster.name] = cluster
        # Verbose output
        if verbose:
            print('Loaded {:d} clusters'.format(len(clusters)), end=' ')
            print('from {:s}'.fromat(clusters_path))
        # Return clusters dictionary
        return clusters

    def make_consensus(self, clusters_path, clusters_acc, verbose=False):
        """ Generate consensus FASTA file

        TODO First, loop through each HMM in mgnifam, create HMM, run hmmemit
        TODO Then, append fasta sequences retireved from hmmemit to consensus
        fasta file

        Args
        clusters_acc (set)          Set of cluster accession numbers
        clusters_path (str)         Path to clusters root directory
        """
        # Initialize timers
        time_beg, time_end = time(), 0.0
        # Initialize consensus file path
        consensus_path = os.path.join(clusters_path, 'CONSENSUS.fa')
        # Open consensus file in write mode
        consensus_file = open(consensus_path, 'w')
        # Loop through each cluster accession number in MGnifam
        for accession in clusters_acc:
            # Make new temporary HMM model file
            hmm_path = NamedTemporaryFile(delete=False, suffix='.hmm').name
            # Retrieve HMM for current cluster accession number
            self.mgnifam_db.make_hmm(accession=accession, path=hmm_path)

            # Define temporary output file
            out_path = NamedTemporaryFile(delete=False, suffix='.hmm').name
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
            self.hmmpress(hmm_path=hmm_path)

            # Get current directory
            hmm_dir = os.path.dirname(hmm_path)
            # Initialize generated file paths
            h3p_path = next(iglob(hmm_dir + '/*.h3p'))
            h3m_path = next(iglob(hmm_dir + '/*.h3m'))
            h3i_path = next(iglob(hmm_dir + '/*.h3i'))
            h3f_path = next(iglob(hmm_dir + '/*.h3f'))
            # Return generated files path
            return h3p_path, h3m_path, h3i_path, h3f_path

        except CalledProcessError:
            raise CalledProcessError('Error: could not run hmmpress scirpt')

    def make_scan(self, *args, **kwargs):
        """ Scan MGnifam

        Scan all MGnifam clusters in root directory, using previously
        generated consensus sequence
        """
        raise NotImplementedError
