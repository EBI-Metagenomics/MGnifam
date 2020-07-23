# Dependencies
from src.msa.transform import Compose, OccupancyTrim, OccupancyFilter
from src.hmm.hmmer import HMMBuild, HMMSearch
from src.dataset import LinClust
from src.dataset import Fasta
from src.disorder import MobiDbLite
from src.disorder import compute_threshold, compute_comp_bias
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
    def run(self, cluster_names, clusters_path, min_jobs=1, max_jobs=100, verbose=False):
        """ Run BATCH pipeline

        Args
        cluster_names (set)         Set of cluster names whose SEED alignment
                                    must be generated
        clusters_path (str)         Path where to store results, an attempt of
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
        # Request jobs
        cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Fill cluster members dictionary
        time_run, cluster_members = benchmark(
            fn=self.search_cluster_members,
            cluster_names=cluster_names,
            client=client,
            verbose=verbose
        )

        # Initialize set of retrieved sequence accession numbers
        sequences_acc = set()
        # Loop through each cluster
        for cluster in cluster_members.values():
            # Update sequence accession numbers
            sequences_acc |= set(cluster)

        # Get number of clusters
        num_clusters = len(cluster_members)
        # Get number of sequences
        num_sequences = sum([len(c) for c in cluster_members.values()])

        # Verbose log
        if verbose:
            print('Retrieved {:d} accession numbers'.format(num_sequences), end=' ')
            print('for {:d} clusters'.format(num_clusters), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Initialize mgnifam dict(sequence acc: fasta entry)
        mgnifam_seq = dict()
        # TODO Fill mgnifam dictionary
        time_run, fasta_sequences = benchmark(
            fn=self.search_fasta_sequences,
            sequences_acc=sequences_acc,
            fasta_dataset=self.mgnifam,
            client=client,
            verbose=verbose
        )

        # Define number of sequences
        num_sequences = len(fasta_sequences)

        # Verbose log
        if verbose:
            print('Retrieved {:d} FASTA'.format(num_sequences), end=' ')
            print('sequences in {:.0f} seconds'.format(time_run))

        # Check compositional bias
        time_run, (bias_below, bias_above) = benchmark(
            fn=self.check_comp_bias,
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print('Compositional bias computed', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

        # Make SEED alignments
        time_run, _ = benchmark(
            fn=self.make_seed_alignments,
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            clusters_path=clusters_path,
            client=client,
            verbose=verbose
        )

        # Verbose log
        if verbose:
            print('SEED sequence alignments generated', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('in {:.0f} seconds'.format(time_run))

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
    def search_cluster_members(self, cluster_names, client, verbose=False):
        """ Search for cluster member sequences

        Args
        cluster_names (set)         Set of cluster names to search against
                                    LinClust dataset
        client (Client)             Dask Client used to schedule jobs
        verbose (bool)              Whether to print verbose log or not

        Return
        (dict)                      Mapping cluster names (keys) to set of
                                    sequence accession numbers (values)
        """
        # Verbose log
        if verbose:
            print('Searching member sequences', end=' ')
            print('for {:d} clusters'.format(len(cluster_names)))

        # Initialize futures
        futures = list()
        # Go through each LinClust dataset chunk
        for chunk in self.linclust:
            # Submit search function to client
            future = client.submit(chunk.search, cluster_names=cluster_names)
            # Store future
            futures.append(future)

        # Retrieve distributed computation results
        results = client.gather(futures)

        # Initialize cluster members dictionary (output)
        cluster_members = dict()
        # Loop through resulting clusters dictionaries
        for i in range(len(results)):
            # Loop through every cluster name in result dictionary
            for cluster_name in results[i].keys():
                # Initialize entry for current cluster
                cluster_members.setdefault(cluster_name, set())
                # Update cluster members list
                cluster_members[cluster_name] |= set(results[i][cluster_name])

        # Verbose log
        if verbose:
            print('Retrieved sequences', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)))

        # Return cluster members
        return cluster_members

    # Search for fasta sequences
    def search_fasta_sequences(self, sequences_acc, fasta_dataset, client, verbose=False):
        """ Retrieve fasra sequences from a Fasta dataset

        Args
        sequences_acc (set)         Set of sequence accession numbers to search
                                    for against given target FASTA dataset
        fasta_dataset (Dataset)     Target FASTA dataset instance, whose search
                                    method will be used in order to find
                                    FASTA entries associated to given
                                    accession numbers
        client (Client)             Dask Client used to schedule jobs
        verbose (bool)              Whether to print verbose log or not

        Return
        (dict)                      Dictionary mapping sequence accession
                                    numbers (keys) to list of fasta entries
                                    (values)
        (int)                       Total length of input FASTA dataset
        """
        # Verbose log
        if verbose:
            print('Searching {:d}'.format(len(sequences_acc)), end=' ')
            print('FASTA sequences')

        # Initialize futures
        futures = list()
        # Loop through each chunk in given FASTA dataset
        for chunk in fasta_dataset:
            # Submit search function against FASTA dataset
            future = client.submit(chunk.search, sequences_acc, ret_length=True)
            # Save future
            futures.append(future)

        # Gather results
        results = client.gather(futures)

        # Initialize FASTA sequences dict(sequence accession: fasta entry)
        fasta_sequences = dict()
        # Initialize length of FASTA file
        fasta_length = 0
        # Go through each result
        for i in range(len(results)):
            # Retrieve either chunk sequences dict and chunk length
            chunk_sequences, chunk_length = results[i]
            # Update cluster sequences dictionary
            fasta_sequences = {**fasta_sequences, **chunk_sequences}
            # Update length of FASTA file
            fasta_length = fasta_length + chunk_length

        # Verbose log
        if verbose:
            print('Retrieved {:d} fasta'.format(len(fasta_sequences)), end=' ')
            print('entries from dataset of size'.format(fasta_length))

        # Return either the sequences dictionary and the dataset length
        return fasta_sequences, fasta_length

    # Check compositional bias
    def check_comp_bias(self, cluster_members, fasta_sequences, client, verbose=False):
        """ Check compositional bias

        Given a cluster members dictionary and a fasta sequences dictionary,
        makes clusters dict(cluster members: fasta sequences). Then makes
        disorder predictions over the sequences: if rate of disordered residues
        is too high, cluster is discarded

        Args
        cluster_members (dict)      Dictionary mappung cluster names to cluster
                                    seqeunces accession numbers
        fasta_sequences (dict)      Dictionary mapping seqeunces accession
                                    numbers to fasta entries
        client (Client)             Dask client used to submit jobs
        verbose (bool)              Wether to print verbose output or not

        Return
        (dict)                      Input cluster entries above compositional
                                    bias threshold
        (dict)                      Input cluster entries below compositional
                                    bias threshold
        """
        # Verbose log
        if verbose:
            print('Computing compositional bias', end=' ')
            print('for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('and {:d} sequences'.format(len(fasta_seqeunces)))

        # Initialize futures dict(cluster name: compositional bias)
        futures = dict()
        # Loop through each cluster
        for cluster_name in cluster_members:

            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = dict()
            # Loop through each sequence accession in cluster
            for acc in cluster_members.get(cluster_name):
                # Save fasta sequence
                cluster_sequences[acc] = fasta_sequences[acc]

            # Run distributed compositional bias computation
            future = client.submit(
                func=self.compute_comp_bias,
                fasta_sequences=cluster_sequences,
                mobidb=self.mobidb,
                pred_threshold=1,  # Set threshold for disorder prediction
                pred_inclusive=True  # Set thresold inclusive
            )

            # Save future
            futures.append(future)

        # Retrieve results
        results = client.gather(futures.values())

        # Initialize set of clusters below and above threshold
        bias_below, bias_above = set(), set()
        # Loop through each cluster name in results
        for i, cluster_name in enumerate(futures.keys()):
            # Get current compositional bias value
            curr_bias = results[i]
            # Case compositional bias threshold is inclusive
            if self.comp_bias_inclusive:
                # Case compositional bias is above threshold
                if curr_bias >= self.comp_bias_threshold:
                    # Add cluster name to the above threshold set
                    bias_above.add(cluster_name)
            # Case compositional bias threshold is exclusive
            elif curr_bias > self.comp_bias_threshold:
                # Add cluster name to the above threshold set
                bias_above.add(cluster_name)
            # Otherwise
            else:
                # Add cluster name to the below threshold set
                bias_below.add(cluster_name)

        # Verbose log
        if verbose:
            print('Compositional bias computed for {:d} clusters'.format(len(cluster_members)), end=' ')
            print('among which {:d} were below'.format(len(bias_below)), end=' ')
            print('threshold of {:f},'.format(self.comp_bias_threshold), end=' ')
            print('while {:d} were above it'.format(len(bias_above)), end=' ')

        # Return both below and above threshold sets
        return bias_below, bias_above

    @staticmethod
    def compute_comp_bias(fasta_sequences, mobidb, pred_threshold=1, pred_inclusive=True):
        """ Compute compositional bias for a cluster

        Given a sequences dict(seqeunce acc: fasta entry), run MobiDbLite
        predictor against them and produces a disordered prediction for each
        residue. A residue is considered to be disordere if it has at least
        <pred_threshold> disordered predictions.

        Once disordered predictiona are computed, compositional bias is
        calculated as rate of predicted disordered residues over number of
        total residues.

        Args
        seqeunces (dict)        Dictionary mapping sequence accession to actual
                                FASTA entry
        mobidb (MobiDbLite)     Instance of MobiDB Lite interface
        pred_threshold (int)    Number of times a residue must be predicted
                                disordered to be considered disordered.
        pred_inclusive (bool)   Wether the prediction threshold should be
                                inclusive or not

        Return
        (float)                 Compositional bias score
        """
        # Make predictions running MobiDB Lite
        predictions = mobidb.run(sequences=fasta_sequences)
        # Apply threshold over predictions
        predictions = compute_threshold(
            sequences=predictions,
            threshold=pred_threshold,
            inclusive=pred_inclusive
        )
        # Compute and return compositional bias
        return compute_comp_bias(predictions)

    # Make SEED alignments
    def make_seed_alignments(self, cluster_members, fasta_sequences, clusters_path, client, verbose=False):
        """ Make SEED alignments

        Given a cluster memebers dictionary associating sequences accession
        numbers to a cluster name and a fasta sequences dictionary associating
        to each sequence accession number its fasta entry, SEED alignment is
        produced through Muscle algorithm and stored at given cluster directory.

        Args
        cluster_members (dict)      Dictionary associating cluster names (keys)
                                    to its member sequences accession numbers
                                    (values)
        fasta_seqeunces (dict)      Dictionary associating sequence accession
                                    numbers (keys) to fasta entries (values)
        clusters_path (str)         Path where cluster folders can be found
        client (str)                Client used to schedule jobs
        verbose (bool)              Wether to show verbose log or not
        """
        # Verbose log
        if verbose:
            print('Making SEED alignments', end=' ')
            print('for {:d}'.format(len(cluster_members)), end=' ')
            print('at {:s}'.format(clusters_path))

        # Initialize futures
        futures = list()
        # Loop through each cluster
        for cluster_name in cluster_members:

            # Define current cluster path
            cluster_path = os.path.join(clusters_path, cluster_name)
            # Ensure that cluster directory exists
            if not os.path.exists(cluster_path):
                # Try to make one
                os.mkdir(cluster_path)

            # Define raw alignment path
            alignment_path = os.path.join(cluster_path, 'SEED_raw')

            # Initialize cluster sequences
            cluster_sequences = list()
            # Loop through each sequence accession in cluster
            for sequence_acc in cluster_members.get(cluster_name):
                # Save fasta sequence
                cluster_sequences.append(fasta_sequences[sequence_acc])

            # Make SEED alignment in cluster directory
            future = client.submit(
                func=self.make_muscle_alignment,
                fasta_sequences=cluster_sequences,
                alignment_path=alignment_path,
                muscle=self.muscle,
            )

            # Store future
            futures.append(future)

        # Retrieve results
        client.gather(futures)

    @staticmethod
    def make_muscle_alignment(fasta_sequences, alignment_path, muscle):
        """ Use MUSCLE to make a multiple sequence alignment

        Args
        fasta_sequences (list)      List of fasta sequences
        alignment_path (str)        Path where alignment will be moved once
                                    done
        muscle (Muscle)             Instance of Muscle interface object
        """
        # Make multiple sequence alignment
        aln = muscle.run(
            sequences=sequences,
            acc_regex='>(\S+)',
            verbose=False
        )
        # Store MSA to disk
        aln.to_aln(out_path=alignment_path)
