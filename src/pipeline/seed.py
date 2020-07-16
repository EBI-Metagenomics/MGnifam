# Dependencies
import numpy as np
import shutil
import time
import sys
import os

from src.utils import get_paths
from src.pipeline.pipeline import Pipeline, Log
from src.dataset import LinClust, Fasta
from src.msa import MSA, Muscle
from src.hmm import HMM, HMMBuild, HMMSearch
from src.disorder import MobiDBLite
import src.disorder as disorder

# Plotting dependencies
import matplotlib.pyplot as plt
import matplotlib
# Do not use interactive plotting
matplotlib.use('Agg')


class SeedPipeline(Pipeline):
    """Seed alignment block
    This block handles all operations needed to produce a SEED alignment for
    each cluster. Then, automatically trims and filters out SEED alignments.
    Finally, discards clusters whose HMM have at least a match into UniProt.
    """

    # Constructor
    def __init__(
        self, linclust_path, mgnify_path, dask_client,
        mobidb_env=os.environ.copy(), mobidb_cmd=['mobidb_lite.py'],
        muscle_env=os.environ.copy(), muscle_cmd=['muscle'],
        hmm_build_env=os.environ.copy(), hmm_build_cmd=['hmmbuild'],
        hmm_search_env=os.environ.copy(), hmm_search_cmd=['hmmsearch']
    ):
        # Retrieve clusters dataset
        self.ds_linclust = LinClust.from_list(get_paths(linclust_path))
        # Retrieve MGnify dataset
        self.ds_mgnify = Fasta.from_list(get_paths(mgnify_path))
        # Initialize a new MobiDB Lite instance
        self.mobidb = MobiDBLite(cmd=mobidb_cmd, env=mobidb_env)
        # Initialize a new Muscle instance
        self.muscle = Muscle(cmd=muscle_cmd, env=muscle_env)
        # Initialize a new hmmbuild instance
        self.hmm_build = HMMBuild()
        # Initialize a new hmmsearch instance
        self.hmm_search = HMMSearch()
        # Set client (Client(LocalCluster) by default)
        self.dask_client = self.get_client()

    # Run Seed alignment (for a single cluster)
    def run(
        self, cluster_names, clusters_dir,
        comp_bias_threshold=0.8, comp_bias_inclusive=True,
        verbose=False
    ):
        """Run SEED and HMM creation block
        Takes as input the <cluster_names> list and retrieves cluster members.
        From cluster members, sequences are retrieved and aligned into a SEED
        alignment. Generally, clusters are assigned a specific folder inside a
        given root folder, namely <clusters_dir>. However, clusters whose
        compositional bias exceed a given threshold are moved from that root
        directory to a BIAS sub directory and excluded from ongoing pipeline.
        Lastly, HMMs are built for SEED alignments in each remaining cluster and
        searched against UniProt dataset: if a match is found, cluster is moved
        to UniProt folder and excluded from ongoing pipeline steps.

        Args
        cluster_names (list)        List of cluster names for which SEED alignment
                                    must be done
        clusters_dir (str)          Root folder for clusters folders and logs
        comp_bias_threshold (float) Threshold for maximum compositional bias a
                                    cluster can have to not be discarded
        verbose (bool)              Wether to print progress to output shell
        """

        # Case output directory does not exits
        if not os.path.exists(clusters_dir):
            # Make directory, this can raise exception
            os.mkdir(clusters_dir)
        # Define log path (only if required)
        log_path = os.path.join(clusters_dir, 'log.json')
        # Initialize new log
        log = Log(log_path=log_path, log_dict={
            # Path to clusters directory
            'clusters_dir': clusters_dir,
            # Time required to run the pipeline entirely
            'tot_time': 0.0,
            # Initial clusters information
            'init_clusters': {
                'num_clusters': 0,
                'num_sequences': 0,
                'avg_sequences': 0.0,
                'avg_bias': 0.0,
                'bias_plot': ''
            },
            # BIAS clusters information
            'bias_clusters': {
                'num_clusters': 0,
                # 'num_sequences': 0,
                'avg_sequences': 0.0,
                'avg_bias': 0.0
            },
            # # UniProt clusters information
            # 'uniprot_clusters': {
            #     'num_clusters': 0,
            #     'num_sequences': 0,
            #     'avg_sequences': 0.0,
            #     'avg_bias': 0.0
            # },
            # Resulting clusters information
            'result_clusters': {
                'num_clusters': 0,
                # 'num_sequences': 0,
                'avg_sequences': 0.0,
                'avg_bias': 0.0
            }
        })

        # Initialize timers
        time_beg, time_end = time.time(), None

        # Retrieve cluster members from clusters .tsv dataset
        cluster_members = self.get_cluster_members(
            cluster_names=cluster_names,
            verbose=verbose
        )
        # Check cluster members (raises error)
        self.check_cluster_members(cluster_members)

        # Define sequences accession numbers list
        sequences_acc = set([
            # Get i-th sequence accession in cluster named n
            cluster_members[n][i]
            # For each cluster in cluster members
            for n in cluster_members
            # For each component of that cluster
            for i in range(len(cluster_members[n]))
        ])

        # Retrieve fasta sequences from MGnify .fa dataset
        fasta_sequences, mgnify_length = self.get_fasta_sequences(
            sequences_acc=sequences_acc,
            verbose=verbose
        )

        # Check for mismatch between cluster sizes and number of inner sequences
        self.check_fasta_sequences(cluster_members, fasta_sequences)

        # Create clusters directories
        self.make_cluster_dirs(cluster_names, clusters_dir)

        # Compute compositional bias
        comp_bias = self.compute_comp_bias(
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            verbose=verbose
        )

        # Define compositional bias plot path
        bias_plot_path = os.path.join(clusters_dir, 'init_comp_bias.png')
        # Plot compositional bias
        self.plot_comp_bias(
            # Compositional bias values
            comp_bias=[*comp_bias.values()],
            # Where to store plot
            out_path=bias_plot_path
        )

        # Update log with initial information
        log({
            # Update initial cluster information
            'init_clusters': {**log.init_clusters, **{
                'num_clusters': len(cluster_names),
                'num_sequences': len(fasta_sequences),
                'avg_sequences': np.mean([
                    len(v) for k, v in cluster_members.items()
                ]),
                'avg_bias': np.mean([
                    v for k, v in comp_bias.items()
                ]),
                'comp_bias_plot': bias_plot_path
            }}
        })

        # Apply compositional bias threshold
        cb_below, cb_above = self.threshold_comp_bias(
            comp_bias=comp_bias,
            threshold=comp_bias_threshold,
            inclusive=comp_bias_inclusive
        )
        # Define bias folder path
        bias_dir = os.path.join(clusters_dir, 'BIAS')
        # Make BIAS folder
        os.mkdir(bias_dir)
        # Move every cluster with compositional bias too high to BIAS folder
        for cluster_name in cb_above:
            # Define current path directory
            cluster_path = os.path.join(clusters_dir, cluster_name)
            # Cluster path in BIAS folder
            cluster_bias_path = os.path.join(bias_dir, cluster_name)
            # Move the cluster folder to BIAS folder
            shutil.move(cluster_path, cluster_bias_path)

        # Update log with BIAS information
        log({
            'bias_clusters': {**log.bias_clusters, **{
                'num_clusters': len(cb_above),
                'avg_sequences': np.mean([
                    len(cluster_members[k]) for k in cb_above
                ]),
                'avg_bias': np.mean([
                    cb_above[k] for k in cb_above
                ])
            }}
        })

        # Discard clusters above the compositional bias threshold
        cluster_members = {
            cluster_name: sequences_acc
            for cluster_name, sequences_acc
            in cluster_members.items()
            if cluster_name in set(cb_below.keys())
        }

        # Make seed alignments
        self.make_sequences_aln(
            cluster_members=cluster_members,
            fasta_sequences=fasta_sequences,
            clusters_dir=clusters_dir,
            verbose=verbose
        )

        # # Make HMM
        # self.build_hmms(
        #     cluster_members=cluster_members,
        #     cluster_dir=cluster_dir,
        #     verbose=verbose
        # )
        #
        # # Search HMM against MGnify
        # self.search_hmms(
        #     cluster_members=cluster_members,
        #     cluster_dir=cluster_dir,
        #     ds_size=mgnify_length,
        #     target_ds=self.ds_mgnify,
        #     verbose=verbose
        # )

        # Update timers
        time_end = time.time()
        # Update log
        log({
            'result_clusters': {**log.result_clusters, **{
                'num_clusters': len(cluster_members),
                'avg_sequences': np.mean([
                    len(cluster_members[k]) for k in cluster_members
                ]),
                'avg_bias': np.mean([
                    v for k, v in cb_above.items()
                ])
            }},
            'tot_time': round(time_end - time_beg, 2)
        })

    # Make clusters directories
    def make_cluster_dirs(self, cluster_names, clusters_dir):
        """Create empty clusters directories
        For each cluster name, creates an associated directory.

        Args
        cluster_names (list)    List of cluster names
        clusters_dir (str)      Path wether to make clusters subdirectories

        Raise
        (FileExistsError)       Could raise an error if file already exists
        (FileNotFoundError)     Could raise an error if permissions are not
                                sufficient or given path does not exist
        """
        # Loop through each cluster name to make clusters directories
        for cluster_name in cluster_names:
            # Define cluster subdir
            cluster_dir = os.path.join(clusters_dir, cluster_name)
            # Check if a subdirectory for current cluster exists
            if not os.path.exists(cluster_dir):
                # Make sub directory
                os.mkdir(cluster_dir)

    # Retrieve cluster members for a batch of clusters
    def get_cluster_members(self, cluster_names, verbose=False):
        # Initialize output dict(cluster name: sequence acc)
        cluster_members = dict()
        # Intialize timers
        time_beg, time_end = time.time(), None

        # Verbose log
        if verbose:
            # Get a summary of cluster names
            print('Retrieving {:d} clusters:\n{}'.format(
                len(cluster_names),  # Number of clusters
                ', '.join(cluster_names)  # Verbose list of cluster names
            ))

        # Submit jobs for searching cluster members
        futures = self.dask_client.map(
            # Search for cluster members in current chunk
            lambda ds: ds.search(cluster_names),
            # List of cluster datasets
            self.ds_linclust
        )
        # Retrieve found clusters (list of lists)
        results = self.dask_client.gather(futures)
        # Loop through resulting clusters dictionaries
        for i in range(len(results)):
            # Loop through every cluster name in result dictionary
            for cluster_name in results[i].keys():
                # Initialize entry for current cluster
                cluster_members.setdefault(cluster_name, [])
                # Update cluster members list
                cluster_members[cluster_name] += results[i][cluster_name]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Get execution time summary
            print('Took {:.0f} seconds to search members of {:d} clusters'.format(
                time_end - time_beg,  # Time required to distributed search
                len(cluster_members.keys())  # Number of clusters
            ))

        # Return cluster members dict
        return cluster_members

    # Check that all clusters have at least one sequence accession
    def check_cluster_members(self, cluster_members):
        # Check for clusters with no members (eventually raise error)
        empty_clusters = [n for n, m in cluster_members.items() if not len(m)]
        # Check if at least one cluster has no member
        if empty_clusters:
            # Raise new error
            raise KeyError('No sequences found for {:d} clusters: {:s}'.format(
                len(empty_clusters),  # Number of clusters with no members
                ', '.join(empty_clusters)  # Verbose string of empty clusters
            ))

    # Retrieve cluster sequences for given sequences accession numbers
    def get_fasta_sequences(self, sequences_acc, verbose=False):
        # Initialize output dict(sequence acc: fasta entry) and total number of searched sequences
        fasta_sequences, fasta_length = dict(), 0
        # Initialize timers
        time_beg, time_end = time.time(), None

        # Submit jobs for searching cluster members sequences
        futures = self.dask_client.map(
            # Search for cluster members in current chunk
            lambda ds: ds.search(sequences_acc, ret_length=True),
            # List of cluster datasets
            self.ds_mgnify
        )
        # Retrieve found sequences (list of dict)
        results = self.dask_client.gather(futures)
        # Loop through resulting sequences dictionaries
        for i in range(len(results)):
            # Get retrieved sequences dict and dataset length
            ret_sequences, ret_length = results[i]
            # Loop through retrieved sequence accession numbers
            for acc in ret_sequences:
                # Update cluster sequences dictionary
                fasta_sequences[acc] = ret_sequences[acc]
            # Update total number of searched sequences
            fasta_length += ret_length

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Log time taken for fasta sequences distributed search
            print('Took {:.0f} seconds to search for {:d} sequences in {:d} sequences'.format(
                time_end - time_beg,  # Time required to distributed search
                len(fasta_sequences),  # Number of retrieved sequences
                fasta_length  # Total number of fasta sequences searched
            ))

        # Return fasta sequences
        return fasta_sequences, fasta_length

    # Check that all the sequences have been retrieved
    def check_fasta_sequences(self, cluster_members, fasta_sequences):
        # Initialize list of clusters with mismatching sequences
        mismatch_clusters = list()
        # Get sequence accession numbers in retrieved fasta sequences
        sequences_acc = set(fasta_sequences.keys())
        # Fill mismatching clusters
        for cluster_name in cluster_members.keys():
            # Define current cluster's seuqence accessionnumbers
            cluster_acc = set(cluster_members[cluster_name])
            # Case some sequence has not been retrieved
            if cluster_acc - sequences_acc:
                # Add cluster name to list of mismatching ones
                mismatch_clusters.append(cluster_name)
        # Case there is at least one mismatching cluster
        if mismatch_clusters:
            # Raise new error
            raise KeyError(
                'Could not retrieve all sequences for {:d} clusters: {}'.format(
                    len(mismatch_clusters),  # Number of mismatching clusters
                    ', '.join(mismatch_clusters)  # Mismatching clusters
                )
            )

    @staticmethod
    def compute_comp_bias_(sequences, mobidb, threshold=1, inclusive=True):
        """Compute compositional bias of a single cluster
        First, calls MobiDB Lite to predict disordered regions
        Then, applies a threshold over minimum number of predictions required
        to define a reidue actually disordered.
        Lastly, computes compositional bias over the thesholded predictions.

        Args
        sequences (dict(str: str))      Dictionary linking each cluster member
                                        sequence accession number to the
                                        associated fasta sequence
        mobidb (disorder.MobiDBLite)    Instance of MobiDB Lite
        threshold (int/str)             Threshold on sum of disorder predictions
                                        in order to define a single residue
                                        actually disordered
        inclusive (bool)                Wether the threshold must be
                                        inclusive (i.e. #pred >= threshold) or
                                        exclusive (i.e. #pred > threshold)

        Return
        (float)                         Compositional bias
        """
        # Make disordered regions predictions
        pred = mobidb.run(sequences=sequences)
        # Apply threshold
        pred = disorder.compute_threshold(
            sequences=pred,
            threshold=threshold,
            inclusive=inclusive
        )
        # Compute compositional bias
        return disorder.compute_comp_bias(pred)

    # Compute compositional bias for every sequence
    def compute_comp_bias(self, cluster_members, fasta_sequences, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None

        # Initialize output dict(cluster name: compositional bias)
        comp_bias = dict()

        # Get cluster names
        cluster_names = [*cluster_members.keys()]

        # Initialize futures list
        futures = list()
        # Loop through every cluster name
        for cluster_name in cluster_names:
            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = dict()
            # Loop through each sequence accession in cluster
            for acc in cluster_members[cluster_name]:
                # Save fasta sequence
                cluster_sequences[acc] = fasta_sequences[acc]
            # Run distributed compositional bias
            futures.append(self.dask_client.submit(
                self.compute_comp_bias_,
                sequences=cluster_sequences,
                mobidb=self.mobidb,
                threshold=1,
                inclusive=True
            ))
        # Get results
        results = self.dask_client.gather(futures)
        # Associate cluster name to each compositional bias
        for i in range(len(cluster_names)):
            # Save compositional bias for current cluster
            comp_bias[cluster_names[i]] = results[i]

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            # Get execution time summary
            print('Took {:.0f} seconds to compute disordered regions for {:d} fasta sequences'.format(
                time_end - time_beg,  # Time required to distributed search
                len(comp_bias)  # Number of fasta sequences
            ))

        # Return dict(cluster name: compositional bias)
        return comp_bias

    def threshold_comp_bias(self, comp_bias, threshold=0.8, inclusive=True):
        """Applies a threshold over compositional bias

        Args
        comp_bias (dict)    Dictionary containing cluster names as keys and
                            compositional bias as values
        threshold (float)   A threshold to use over compositional bias to
                            eventually discard some clusters

        Return
        (dict)              Same dict as input, but containing only entries
                            below given threshold
        (dict)              Same dict as input, but containing only entries
                            above given threshold
        """
        # Intialize output dictionaries
        below_threshold, above_threshold = dict(), dict()
        # loop through all entries in input dict
        for k, v in comp_bias.items():
            # Check if compositional bias value exceeds threshold
            exceed = (v >= threshold) if inclusive else (v > threshold)
            # Case compositional bias value exceeds threshold
            if exceed:
                # Assign current cluster to above threshold ones
                above_threshold[k] = v
            # Otherwise
            else:
                # Assign current cluster to below threshold ones
                below_threshold[k] = v
        # Return output dictionaries
        return below_threshold, above_threshold

    # Make a plot of compositional bias
    def plot_comp_bias(self, comp_bias, out_path):
        # Plot compositional bias distribution
        fig, ax = plt.subplots(figsize=(10, 5))
        # Make plot
        ax.set_title('Compositional Bias distribution')
        ax.set_xlabel('Compositional Bias')
        ax.set_ylabel('Number of clusters')
        ax.hist(comp_bias, density=False, bins=100)
        ax.set_xlim(left=0.0, right=1.0)
        # Save plot
        plt.savefig(out_path)
        # Close plot
        plt.close()

    @staticmethod
    def make_sequences_aln_(sequences, muscle, aln_path):
        """Aligns a list of fasta sequences
        Takes as input a list of sequences and aligns them, then store the
        alignment as a .aln file.

        Args
        sequences (dict(str: str))
        muscle (msa.Muscle)         Multiple Sequence Alignment (MSA) algorithm
                                    used for generating output .aln file
        path (str)                  Path to output multiple sequence alignment
                                    file
        """
        # Make multiple sequence alignment
        aln = muscle.run(sequences=sequences, acc_regex='>(\S+)', verbose=False)
        # Store MSA to disk
        aln.to_aln(out_path=aln_path)
        # Free memory
        del aln

    # Make multiple sequence alignment
    def make_sequences_aln(self, cluster_members, fasta_sequences, clusters_dir, verbose=False):

        # Initialize timers
        time_beg, time_end = time.time(), None
        # Verbose log
        if verbose:
            print('Making multiple sequence alignments for {:d} clusters...'.format(
                len(cluster_members)
            ))

        # Get cluster names
        cluster_names = [*cluster_members.keys()]

        # Initialize futures list
        futures = list()
        # Loop through every cluster name
        for cluster_name in cluster_names:
            # Define a cluster sequences dict(sequence acc: fasta entry)
            cluster_sequences = list()
            # Loop through each sequence accession in cluster
            for acc in cluster_members[cluster_name]:
                # Save fasta sequence
                cluster_sequences.append(fasta_sequences[acc])
            # Define path to SEED alignment
            seed_path = os.path.join(clusters_dir, cluster_name, 'SEED')
            # Run distributed compositional bias
            futures.append(self.dask_client.submit(
                self.make_sequences_aln_,
                sequences=cluster_sequences,
                muscle=self.muscle,
                aln_path=seed_path
            ))
        # Get results
        self.dask_client.gather(futures)

        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            print('Took {:.0f} seconds to make {:d} multiple sequence alignments'.format(
                time_end - time_beg,
                len(cluster_members)
            ))

    # # Create Hidden Markov Models (HMMs)
    # def build_hmms(self, cluster_members, cluster_dir, verbose=False):
    #
    #     # Initialize timers
    #     time_beg, time_end = time.time(), None
    #     # Verbose log
    #     if verbose:
    #         print('Building Hidden Markov Models (HMMs) for {:d} clusters...'.format(
    #             len(cluster_members)
    #         ))
    #
    #     # Intialize futures container
    #     futures = list()
    #     # Loop through each cluster name
    #     for cluster_name in cluster_members:
    #         # Define path to current cluster
    #         cluster_path = os.path.join(cluster_dir, cluster_name)
    #         # Define path to current cluster seed alignment
    #         seed_path = os.path.join(cluster_path, 'SEED')
    #         # Define path to output Hidden Markov Model (HMM)
    #         hmm_path = os.path.join(cluster_path, 'HMM')
    #         # Run HMM build (distributed)
    #         futures.append(self.dask_client.submit(
    #             # Function to be submitted
    #             self.hmm_build.run,
    #             # Parameters to be feed to submitted function
    #             msa_path=seed_path,  # Set path to input seed alignment
    #             out_path=hmm_path,  # Set path to output file
    #             name=cluster_name  # Set HMM name as current cluster name
    #         ))
    #
    #     # Retrieve results
    #     self.dask_client.gather(futures)
    #
    #     # Update timers
    #     time_end = time.time()
    #     # Verbose log
    #     if verbose:
    #         print('Took {:.0f} seconds to build {:d} Hidden Markov Models (HMMs)'.format(
    #             time_end - time_beg,
    #             len(cluster_members)
    #         ))
    #
    # @staticmethod
    # def search_hmms_(hmm_search, hmm_paths, target_path, num_cpus=1):
    #     """Search an HMM file/library against a target dataset
    #     Target dataset is compressed, so it gets uncompressed firstly.
    #     Then, domtblout is retrieved by running hmmsearch.
    #
    #     Args
    #     hmm_search (HMMSearch)      Already initialized wrapper for hmmsearch
    #                                 script
    #     hmm_paths (list)            List of paths to HMM files to concatenate
    #                                 in a sinfle HMM library
    #     target_path (str)           Path to target dataset, against wich given
    #                                 HMM library will be searched
    #     num_cpus (int)              Number of CPUs to use in hmmsearch run
    #
    #     Return
    #     (list)                      List containing domtblout table rows
    #     """
    #
    #     # Initialize current HMM library empty file
    #     hmm_lib_path = NamedTemporaryFile(suffix='.HMM', delete=False).name
    #     # Fill current HMM library file
    #     HMM.concat(in_paths=hmm_paths, out_path=hmm_lib_path)
    #
    #     # Initialize a temporary input and output files
    #     in_path = NamedTemporaryFile(suffix='.fasta', delete=False).name
    #     out_path = NamedTemporaryFile(delete=False).name
    #
    #     # Uncompress input file (if compressed)
    #     gunzip(in_path=target_path, out_path=in_path)
    #
    #     # Search given hmm aganist given dataset using hmmsearch
    #     hmm_search.run(
    #         # Path to searched HMM model/library file
    #         hmm_path=hmm_lib_path,
    #         # Path to sequences target dataset
    #         ds_path=in_path,
    #         # Path to per-domain output table
    #         domtblout_path=out_path,
    #         # Define number of CPUs allocated
    #         num_cpus=num_cpus
    #     )
    #
    #     # Initialize domtblout
    #     domtblout = list()
    #     # Open output temporary file
    #     with open(out_temp.name) as out_file:
    #         # Save domtblout as list of rows
    #         for row in iter_domtblout(out_file):
    #             # Store current row
    #             domtblout.append(row)
    #
    #     # Remove temporary input and output files
    #     os.remove(in_path)
    #     os.remove(out_path)
    #
    #     # Return per-domain table
    #     return domtblout
    #
    # # Search Hidden Markov Models (HMMs) against a dataset
    # def search_hmms(self, cluster_members, cluster_dir, target_ds, ds_size, verbose=False):
    #     """Search multiple HMMs against given dataset
    #     First, goes through each HMM path and checks if the required resources
    #     fit into allocated ones. If not, skips the cluster, otherwise adds it to
    #     HMM library which will be checked against the whole target dataset.
    #     Finally, HMM library is run against each chunk in the target dataset
    #     in a distributed manner.
    #
    #     Args
    #     cluster_members (iterable)  Contains names of clusters whose HMM must be
    #                                 searched against given dataset
    #     cluster_dir (str)           Root directory of clusters sub-directories
    #     target_ds (list)            List of dataset chunks against which the
    #                                 clusters Hidden Markov Models must be run
    #                                 in a distributed way
    #     ds_size (int)               Number of entries in input dataset, in total
    #     verbose (bool)              Wether to return or not verbose output
    #
    #     Return
    #     (list)                      Concatenated domtblout table, containing
    #                                 all matches between HMMs and sequences
    #     """
    #
    #     # List of HMMs
    #     hmm_lib = list()
    #     # Initialize maximum available resources
    #     max_memory, max_cpus = 8, 8
    #     # Loop through each cluster names
    #     for cluster_name in cluster_members:
    #         # Define full path to current HMM file
    #         hmm_path = os.path.join(cluster_dir, cluster_name, 'HMM')
    #         # Load HMM from file
    #         hmm_obj = HMM.from_file(hmm_path)
    #         # Get HMM length
    #         hmm_len = hmm_obj.length
    #         # Compute required memory resources
    #         num_cpus = HMMSearch.get_resources(
    #             max_memory=max_memory,  # Maximum memory in Gb
    #             model_len=hmm_len,  # Model length
    #             longest_seq=4e04,  # Longest sequence
    #             max_cpus=max_cpus  # Set maximum number of CPUs
    #         )
    #         # Case current HMM does not fit allocated resources
    #         if not num_cpus:
    #             continue  # Skip iteration
    #         # Otherwise, update maximum number of CPUs
    #         max_cpus = min(num_cpus, max_cpus)
    #         # Otherwise, add HMM path to HMM library
    #         hmm_lib.append(hmm_path)
    #
    #     # # Try running HMM search commands
    #     # try:
    #
    #     # Initialize futures container
    #     futures = list()
    #     # Loop through each dataset in <against_ds> list
    #     for i in range(len(target_ds)):
    #         # Distribute HMM search against target dataset
    #         futures.append(self.dask_client.submit(
    #             # Function to submit
    #             self.search_hmms_,
    #             # Function parameters
    #             hmm_search=self.hmm_search,  # HMMSearch object
    #             hmm_path=hmm_lib_path,  # Target HMM library
    #             target_ds=target_ds[i].path,  # Target dataset path
    #             num_cpus=max_cpus  # Number of CPUs allocated
    #         ))
    #
    #     # Retrieve results
    #     results = self.dask_client.gather(futures)
    #
    #     # # Remove input HMM library
    #     # os.remove(hmm_lib_path)
    #
    #     # # Subprocess error
    #     # except CalledProcessError as err:
    #     #     # Debug
    #     #     print('hmmsearch exited with code', err.returncode, file=sys.stderr)
    #     #     print('cmd:', err.cmd, file=sys.stderr)
    #     #     print('stderr:', file=sys.stderr)
    #     #     print(err.stderr, file=sys.stderr)
    #     #
    #     # # Client error
    #     # except Exception as err:
    #     #     # Debug
    #     #     traceback.print_exc()
    #     #
    #     # # After exception
    #     # finally:
    #     #     # Exit with error
    #     #     sys.exit(1)
    #
    #     # Concatenate all results in a single table
    #     domtblout = [
    #         # Get j-th domtblout row of i-th dataset chunk
    #         results[i][j]
    #         # Loop through each dataset chunk result
    #         for i in range(len(results))
    #         # Loop through each row in i-th dataset chunk result
    #         for j in range(len(results[i]))
    #     ]
    #
    #     # Debug
    #     print('Retrieved {:d} rows from domtblout'.format(len(results)))
    #
    # # Wrapper for searching given HMMs against Mgnify dataset
    # def search_hmms_mgnify(self, *args, **kwargs):
    #     # Call generic function
    #     return self.search_hmms(*args, **kwargs, target_ds=self.ds_mgnify)
