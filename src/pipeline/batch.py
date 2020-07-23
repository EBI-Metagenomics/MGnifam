# Dependencies
from src.msa.transform import Compose, OccupancyTrim, OccupancyFilter
from src.hmm.hmmer import HMMBuild, HMMSearch
from src.dataset import LinClust
from src.dataset import Fasta
from src.disorder import MobiDbLite
from src.msa import Muscle
from time import time


class BatchPipeline(Pipeline):
    """ Make a batch of LinCLust clusters

    This pipeline takes as input a list of cluster names, for each cluster
    retrieves its member sequences, makes a SEED alignment through Muscle
    multiple sequence alignment algorithm, trims it automatically (and
    eventually discards it, if not valid)
    """

    # Constructor
    def __init__(
        # Pipeline parameters, required to handle job scheduling
        self, cluster_type, cluster_kwargs,
        # Path to datasets (fixed)
        linclust_path, mgnifam_path, uniprot_path,
        # Compositional bias threshold settings
        comp_bias_threshold=0.2, comp_bias_inclusive=True,
        # Automatic trimming settings
        trim_threshold=0.4, trim_inclusive=True,
        filter_thresold=0.5, filter_inclsive=True,
        # Post trimming settings
        seed_min_width=1, seed_min_height=1,
        # Search against UniProt settings
        uniprot_e_value=0.01, uniprot_z_score=None,
        # Command line arguments
        mobidb_cmd=['mobidb_lite.py'],  # Path to MobiDB Lite predictor
        muscle_cmd=['muscle'],  # Path to Muscle alignment algorithm
        hmm_build_cmd=['hmmbuild'],  # Path to hmmbuild script
        hmm_search_cmd=['hmmsearch'],  # Path to hmmsearch script
        # Environmental variables
        env=os.environ.copy()
    ):
        # Call parent constructor
        super().__init__(cluster_type, cluster_kwargs)
        # Save datasets path
        self.linclust = LinClust.from_str(linclust_path)
        self.mgnifam = Fasta.from_str(mgnifam_path)
        self.uniprot = Fasta.from_str(uniprot_path)
        # Save compositional bias parameters
        self.comp_bias_threshold = comp_bias_threshold
        self.comp_bias_inclusive = comp_bias_inclusive
        # Save automatic trimming settings
        self.trim_threshold = trim_threshold
        sefl.trim_inclusive = trim_inclusive
        self.filter_threshold = filter_threshold
        self.filter_inclusive = filter_inclusive
        # Save post trimming settings
        self.seed_min_width = seed_min_width
        self.seed_min_height = seed_min_height
        # Define automatic trimming
        self.seed_transform = Compose([
            # Exclude regions outside N- and C- terminal
            OccupancyTrim(
                threshold=trim_threshold,
                inclusive=trim_inclusive
            ),
            # Exclude sequences with less than half occupancy
            OccupancyFilter(
                threshold=filter_threshold,
                inclusive=filter_inclusive
            )
        ])
        # Save UniProt hmmsearch settings
        self.uniprot_e_value = uniprot_e_value
        self.uniprot_z_score = uniprot_z_score
        # Define scripts interfaces
        self.mobidb = MobiDbLite(cmd=mobidb_cmd, env=env)
        self.muscle = Muscle(cmd=muscle_cmd, env=env)
        self.hmm_build = HMMBuild(cmd=hmm_build_cmd, env=env)
        self.hmm_search = HMMSearch(cmd=hmm_search_cmd, env=env)
        # Store environmental variables
        self.env = env

    # Run the pipeline
    def run(self, cluster_names, out_path, min_jobs=1, max_jobs=100, verbose=True):
        """ Run BATCH pipeline

        Args
        cluster_names (set)         Set of cluster names whose SEED alignment
                                    must be generated
        out_path (str)              Path where to store results, an attempt of
                                    making the output directory is made if it
                                    does not exist yet.

        Return

        Raise
        """
        # TODO
        raise NotImplementedError

        # Initialize timers
        time_beg, time_end = time(), 0.0

        # Verbose log
        if verbose:
            print('Running BATCH pipeline', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('at {:s}'.format(out_path))

        # Initialize new cluster
        cluster, client = self.get_client({'cores': 1})

        # Initialize cluster members dict(cluster name: sequence acc)
        cluster_members = dict()
        # Fill cluster members dictionary
        time_run, cluster_members = benchmark(
            fn=self.search_cluster_members,
            cluster_names=cluster_names
        )

        # Initialize mgnifam dict(sequence acc: fasta entry)
        mgnifam_seq = dict()
        # TODO Fill mgnifam dictionary

        # TODO Check compositional bias

        # TODO Make SEED alignments

        # TODO Trim SEED alignments

        # TODO Make HMM models

        # TODO Search HMM models against UniProt

        # Update timers
        time_end = time()
        time_tot = time_end - time_beg

        # Verbose log
        if verbose:
            print('BATCH pipeline done', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)), end=' ')
            print('in {:.0f} seconds'.format(time_tot))

        # Return time taken to run and set of cluster names
        return time_tot, cluster_names

    # Search for cluster members
    def search_cluster_members(self, cluster_names, client):
        """ Search for cluster member sequences

        Args
        cluster_names (set)         Set of cluster names to search against
                                    LinClust dataset
        """
