# Dependencies
from tempfile import NamedTemporaryFile
from time import time
import sys
import os

# Custom dependencies
from src.utils import get_paths, as_bytes, is_gzip, gunzip
from src.pipeline.pipeline import Pipeline, Log
from src.hmm import HMM, HMMBuild, HMMSearch
from src.hmm import iter_domtblout
from src.dataset import Fasta


class HMMPipeline(Pipeline):

    def __init__(
        self, cluster_type, cluster_kwargs, env=os.environ.copy(),
        hmm_build_cmd=['hmmbuild'], hmm_search_cmd=['hmmsearch']
    ):
        # Call parent constructor
        super().__init__(cluster_type, cluster_kwargs)
        # Save environmental variables
        self.env = env
        # Initialize a new hmmbuild instance
        self.hmm_build = HMMBuild(cmd=hmm_build_cmd, env=self.env)
        # Initialize a new hmmsearch instance
        self.hmm_search = HMMSearch(cmd=hmm_search_cmd, env=self.env)

    def run(
        self, cluster_names, clusters_dir, target_path, e_value=0.01, z_score=None,
        min_jobs=0, max_jobs=1, log_path='', verbose=False
    ):
        """Handle HMM models

        Goes through clusters in a given directory, builds HMM models and then
        searches each HMM against a given dataset. Discards those matches whose
        E-value exceed a given threshold.

        Args
        cluster_names (list)        List of cluster names, to which a subfolder
                                    is associated
        clusters_dir (str)          Path where cluster directories can be found
        target_path (str/list)      Path (paths if chunked) to target dataset
        e_value (float)             Maximum allowed e-value to consider an
                                    hmmsearch match as valid
        z_score (int)               Dimension of target dataset
        min_jobs, max_jobs (int)    Minumum and maximum nuber of jobs which can
                                    be scheduòed at the same time
        verbose (bool)              Wether to print out verbose logs or not
        """

        # Intialize timer
        tot_time = time()
        # Define new log
        log = Log(log_path=log_path, log_dict={
            'tot_time': 0.0,
            'num_clusters': len(cluster_names),
            'target_path': target_path,
            'e_value': e_value,
            'z_score': z_score,
            'build_time': 0.0,
            'search_time':0.0
        })

        # Instantiate new client for running HMMBuild
        client = self.get_client({'cores': 1, 'memory': '2 GB'})
        # Require some jobs
        client.cluster.adapt(minimum=min_jobs, maximum=max_jobs)

        # Get list of target path to FASTA dataset
        target_paths = get_paths(target_path)
        # Case z score is not set
        if not z_score:
            # Initialize timer
            z_time = time()
            # Get fasta datasets
            fasta = Fasta.from_list(target_paths)
            # Get each chunk size
            futures = client.map(len, fasta)
            # Retrieve results
            results = client.gather(futures)
            # Compute z score
            z_score = sum(results)
            # Update log
            log({
                'z_time': time() - z_time,
                'z_score': z_score
            })

        # Initialize build time
        build_time = time()
        # Make HMMs
        self.build_hmms(
            cluster_names=cluster_names,
            clusters_dir=clusters_dir,
            client=client,
            verbose=verbose
        )
        # Close client
        client.shutdown()
        # Update log
        log({'build_time': round(time() - build_time, 2)})

        # Initialize number of cores
        max_cores = self.cluster_kwargs.get('cores', 1)
        # Initialize memory
        max_memory = as_bytes(self.cluster_kwargs.get('memory', '1 GB'))

        # Initialize HMM library
        hmm_lib_path = list()
        # Find mmaximum usable number of cores
        for cluster_name in cluster_names:

            # Define path to current cluster
            cluster_path = os.path.join(clusters_dir, cluster_name)
            # Define path to HMM
            hmm_path = os.path.join(cluster_path, 'HMM')
            # Load HMM model from file
            hmm_obj = HMM.from_file(hmm_path)
            # Retrieve current HMM model attributes
            hmm_name = hmm_obj.name
            hmm_len = hmm_obj.length

            # Verbose
            if verbose:
                print('HMM model {:s} has length {:d}...'.format(
                    hmm_name,  # Name of the HMM model
                    hmm_len  # Length of the HMM model
                ), flush=False, end='')
            # Get maximum number of cores
            cores = HMMSearch.get_cpus(
                max_memory=max_memory,  # Maximum avaliable memory per job
                max_cpus=max_cores,  # Maximum number of CPUs
                model_len=hmm_len,  # HMM model length
                longest_seq=4e04  # Longest sequence length
            )

            # Case job does not fit memory
            if not cores:
                # Verbose
                if verbose:
                    print('discarded due to requirements')
                # Skip iteration
                continue

            # Verbose
            if verbose:
                print('added to library, with {:d} cores maximum'.format(cores))
            # Update max cores
            max_cores = min(max_cores, cores)
            # Add current model to library
            hmm_lib_path.append(hmm_path)

        # Initialize build time
        search_time = time()
        # Instantiate new client for running
        client = self.get_client({'cores': max_cores})
        # Require some jobs
        client.cluster.adapt(minimum=min_jobs, maximum=max_jobs)
        # Run hmmsearch, retrieve clusters which match sequences in dataset
        matching_clusters = self.search_hmms(
            hmm_paths=hmm_lib_path,
            target_paths=target_paths,
            e_value=e_value,
            z_score=z_score,
            client=client,
            cores=max_cores,
            verbose=verbose
        )
        # Close client
        client.shutdown()
        # Update log
        log({'search_time': round(time() - search_time, 2)})

        # Define path to Uniprot directory
        uniprot_path = os.path.join(clusters_dir, 'Uniprot')
        # Make directory
        os.mkdir(uniprot_path)
        # Loop through every cluster matching target dataset
        for cluster_name in matching_clusters:
            # Define cluster path
            cluster_path = os.path.join(clusters_dir, cluster_name)
            # Case cluster path does not exist
            if not os.path.exists(cluster_path):
                continue  # Skip it
            # Otherwise, move it to Uniprot folder
            shutil.move(cluster_path, os.path.join(uniprot_path, cluster_name))

        # Update log
        log({'tot_time': round(time() - tot_time, 2)})

    # Create Hidden Markov Models (HMMs)
    def build_hmms(self, cluster_names, clusters_dir, client, verbose=False):

        # Initialize timers
        time_beg, time_end = time(), None
        # Verbose log
        if verbose:
            print('Building Hidden Markov Models (HMMs) for {:d} clusters...'.format(
                len(cluster_names)
            ))

        # Intialize futures container
        futures = list()
        # Loop through each cluster name
        for cluster_name in cluster_names:
            # Define path to current cluster
            cluster_path = os.path.join(clusters_dir, cluster_name)
            # Define path to current cluster seed alignment
            seed_path = os.path.join(cluster_path, 'SEED')
            # Define path to output Hidden Markov Model (HMM)
            hmm_path = os.path.join(cluster_path, 'HMM')
            # Run HMM build (distributed)
            futures.append(client.submit(
                # Function to be submitted
                self.hmm_build.run,
                # Parameters to be feed to submitted function
                msa_path=seed_path,  # Set path to input seed alignment
                out_path=hmm_path,  # Set path to output file
                name=cluster_name  # Set HMM name as current cluster name
            ))

        # Retrieve results
        client.gather(futures, 'raise')

        # Update timers
        time_end = time()
        # Verbose log
        if verbose:
            print('Took {:.0f} seconds to build {:d} Hidden Markov Models (HMMs)'.format(
                time_end - time_beg,
                len(cluster_names)
            ))

    @staticmethod
    def search_hmms_(
        hmm_search, hmm_paths, target_path, e_value, z_score, num_cpus=1
    ):
        """Search an HMM file/library against a target dataset
        Target dataset is compressed, so it gets uncompressed firstly.
        Then, domtblout is retrieved by running hmmsearch.

        Args
        hmm_search (HMMSearch)      Already initialized wrapper for hmmsearch
                                    script
        hmm_paths (list)            List of paths to HMM files to concatenate
                                    in a sinfle HMM library
        target_path (str)           Path to target dataset, against wich given
                                    HMM library will be searched
        e_value (float)             E-value threshold defining significant matches
        z_score (int)               Dimension of target dataset
        num_cpus (int)              Number of CPUs to use in hmmsearch run

        Return
        (set)                       Set containing HMM model names which have at
                                    least one match in target dataset
        """

        # Initialize current HMM library empty file
        hmm_lib_path = NamedTemporaryFile(suffix='.HMM', delete=False).name
        # Fill current HMM library file
        HMM.concat(in_paths=hmm_paths, out_path=hmm_lib_path)

        # Case input file is compressed
        if is_gzip(target_path):
            # Unzip target path and overwrite given target dataset
            target_path = gunzip(target_path, out_suffix='.fasta')

        # Initialize domtblout path
        domtblout_path = NamedTemporaryFile(delete=False).name

        # Search given hmm aganist given dataset using hmmsearch
        hmm_search.run(
            # Path to searched HMM model/library file
            hmm_path=hmm_lib_path,
            # Path to sequences target dataset
            target_path=target_path,
            # Path to per-domain output table
            domtblout_path=domtblout_path,
            # Define number of CPUs allocated
            num_cpus=num_cpus,
            # Set e-value
            dom_e=e_value,
            # Set z-score
            dom_z=z_score
        )

        # Initialize domtblout
        hmm_matches = set()
        # Open output temporary file
        with open(domtblout_path) as file:
            # Save domtblout as list of rows
            for row in iter_domtblout(file):
                # Get HMM model name, stored in 3rd column
                hmm_matches.add(row[3][1])
        # Remove temporary output file
        os.remove(domtblout_path)
        # Return HMM model names matching target dataset
        return hmm_matches

    # Search Hidden Markov Models (HMMs) against a dataset
    def search_hmms(
        self, hmm_paths, target_paths, e_value, z_score, client, cores,
        verbose=False
    ):
        """Search multiple HMMs against given dataset
        First, goes through each HMM path and checks if the required resources
        fit into allocated ones. If not, skips the cluster, otherwise adds it to
        HMM library which will be checked against the whole target dataset.
        Finally, HMM library is run against each chunk in the target dataset
        in a distributed manner.

        Args
        hmm_paths (list)            List of HMM paths which musth be checked
                                    against target dataset
        target_path (list)          List of target datasets paths (fasta files)
        e_value (float)             E-value threshold defining significant matches
        z_score (int)               Dimension of target dataset
        client (Client)             Client to schedule parallel hmmsearch jobs
        verbose (bool)              Wether to return or not verbose output

        Return
        (set)                       Set of HMM names (by default named according
                                    to belonging cluster) which have at least a
                                    significan match with target dataset
        """
        # Try to run distributed search
        try:

            # Initialize futures container
            futures = list()
            # Loop through each dataset in <against_ds> list
            for i in range(len(target_paths)):
                # Define current target path
                target_path = target_paths[i]
                # Distribute HMM search against target dataset
                futures.append(client.submit(
                    # Function to submit
                    self.search_hmms_,
                    # Function parameters
                    hmm_search=self.hmm_search,  # HMMSearch object
                    hmm_paths=hmm_paths,  # Target HMM library
                    target_path=target_path,  # Target dataset paths
                    num_cpus=cores,  # Number of CPUs allocated
                    e_value=e_value,  # E-value threshold
                    z_score=z_score  # Dataset size
                ))

            # Retrieve results
            results = client.gather(futures, 'raise')

        # Intercept command line error
        except CalledProcessError as err:
            # Raise new error
            raise Exception(err.stderr)

        # Define set of matching HMMs
        hmm_matches = set()
        # Loop through each computation result
        for i in range(len(results)):
            # Add retrieved set of matches
            hmm_matches = hmm_matches | results[i]
        # Verbose log
        if verbose:
            print('There are {:d} HMMs matching target dataset'.format(
                len(hmm_matches)
            ))
        # Return name of HMM matching
        return hmm_matches

        # # Remove input HMM library
        # os.remove(hmm_lib_path)

        # # Subprocess error
        # except CalledProcessError as err:
        #     # Debug
        #     print('hmmsearch exited with code', err.returncode, file=sys.stderr)
        #     print('cmd:', err.cmd, file=sys.stderr)
        #     print('stderr:', file=sys.stderr)
        #     print(err.stderr, file=sys.stderr)
        #
        # # Client error
        # except Exception as err:
        #     # Debug
        #     traceback.print_exc()
        #
        # # After exception
        # finally:
        #     # Exit with error
        #     sys.exit(1)

        # # Concatenate all results in a single table
        # domtblout = [
        #     # Get j-th domtblout row of i-th dataset chunk
        #     results[i][j]
        #     # Loop through each dataset chunk result
        #     for i in range(len(results))
        #     # Loop through each row in i-th dataset chunk result
        #     for j in range(len(results[i]))
        # ]
        #
