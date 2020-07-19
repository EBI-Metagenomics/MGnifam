# Dependencies
from src.pipeline.pipeline import Pipeline, Log
from src.pipeline.seed import SeedPipeline
from src.pipeline.hmm import HMMPipeline
from src.dataset import Fasta
from time import time
import sys
import os
import re

class ReleasePipeline(Pipeline):
    """Make new MGnifam release

    This pipeline containse either the seed pipeline (mgseed script) and the
    hmm pipeline (pfbuild and pfmake scripts). It does not handle LSF cluster
    itself, but demands it to the inner pipelines insted.
    """

    # Constructor
    def __init__(
        self, cluster_type, cluster_kwargs, linclust_paths, mgnify_paths, uniprot_paths,
        occupancy_trim_threshold=0.4, occupancy_trim_inclusive=True,
        occupancy_filter_threshold=0.5, occupancy_filter_inclusive=True,
        mobidb_cmd=['mobidb_lite.py'], muscle_cmd=['muscle'],
        hmm_build_cmd=['hmmbuild'], hmm_search_cmd=['hmmsearch'],
        hmm_align_cmd=['hmmalign'], env=environ.copy()
    ):
        # Call parent constructor
        super().__init__(cluster_type, cluster_kwargs)
        # Set environmental variables
        self.env = env
        # Set path to datasets (could be chunked)
        self.linclust_paths = linclust_paths
        self.mgnify_paths = mgnify_paths
        self.uniprot_paths = uniprot_paths
        # Set inner SEED pipeline
        self.seed_pipeline = SeedPipeline(
            cluster_type=cluster_type, cluster_kwargs=cluster_kwargs,
            occupancy_trim_threshold=occupancy_trim_threshold,
            occupancy_trim_inclusive=occupancy_trim_inclusive,
            occupancy_filter_threshold=occupancy_filter_threshold,
            occupancy_trim_inclusive=occupancy_filter_inclusive
        )
        # Set inner HMM pipeline
        self.hmm_pipeline = HMMPipeline(
            cluster_type=cluster_type, cluster_kwargs=cluster_kwargs,
            hmm_build_cmd=hmm_build_cmd, hmm_search_cmd=hmm_search_cmd,
            hmm_align_cmd=hmm_align_cmd
        )

    # Execute pipeline
    def run(
        self, clusters_names, release_dir, batch_size=1000, max_clusters=None,
        min_jobs=1, max_jobs=100
    ):
        """Run the pipeline

        Args
        cluster_names (iterable)    Name of clusters which must run through
                                    the pipeline
        realease_dir (str)          Path to directory where clusters will be
                                    stored
        batch_size (int)            Number of clusters to run at once
        max_clusters (int)          Maximum number of clusters to run in total
        min_jobs (int)              Minimum number of jobs to spawn
        max_jobs (int)              Maximum number of jobs to spawn
        """
        # Initialize log
        log = Log(log_path=log_path, log_dict={
            'tot_time': 0.0,  # Time (seconds) it took to run
            'uniprot_len': 0,  # length of uniprot dataset
            'uniprot_time': 0.0,  # Time required to get uniprot length
            'mgnify_len': 0,  # Length of mgnify dataset
        })
        # Initialize beginning timer
        tot_time = time()
        # Initialize release directory
        os.mkdir(release_dir)

        # Initialize timer
        uniprot_time = time()
        # Get UniProt length
        uniprot_len = self.fasta_len(
            target_paths=self.uniprot_paths,
            min_jobs=min_jobs,
            max_jobs=max_jobs
        )
        # Update log
        log({
            'uniprot_time': round(time() - uniprot_time, 2),
            'uniprot_len': uniprot_len
        })

        # Initialize timer
        mgnify_time = time()
        # Get MGnify length
        mgnify_len = self.fasta_len(
            target_path=self.mgnify_paths,
            min_jobs=min_jobs,
            max_jobs=max_jobs
        )
        # Update log
        log({
            'mgnify_time': round(time() - mgnify_time, 2),
            'mgnify_len': mgnify_len
        })

        # Define linclust files iterator
        batch_iter = self.iter_clusters(
            clusters_paths=clusters_paths,  # Path to clusters input file(s)
            batch_size=batch_size,  # Size of a sinclue clusters batch
            max_clusters=max_clusters  # Maximum number of iterable clusters
        )
        # Go through every batch
        for batch_index, cluster_names in batch_iter:
            # Define batch path
            batch_path = os.path.join(release_dir, 'batch_{:d}'.format(batch_index))
            # Make batch directory
            os.mkdir(batch_path)

            # Run seed pipeline, retireve log
            self.mgseed(
                cluster_names=cluster_names,
                clusters_dir=batch_path,
                min_jobs=min_jobs,
                max_jobs=max_jobs
            )

            # Run pfbuild pipeline
            self.pfbuild(
                cluster_names=cluster_names,
                clusters_dir=batch_path,
                target_len=uniprot_len,
                min_jobs=min_jobs,
                max_jobs=max_jobs
            )

            # Run mfbuild pipeline
            self.mfbuild(
                cluster_names=cluster_names,
                clusters_dir=batch_path,
                target_len=uniprot_len,
                min_jobs=min_jobs,
                max_jobs=max_jobs
            )

        # Define build path
        build_path = os.path.join(release_dir, 'build')
        # Make build directory
        os.mkdir(build_path)

        # Define kept clusters
        kept_clusters = os.path.join(release_dir, 'batch_*', 'MGYP*')
        # Go through each kept cluster and copy it into build
        for cluster_path in iglob(kept_clusters):
            # Get cluster name
            cluster_name = os.path.basename(cluster_path)
            # Copy cluster to build path
            shutil.copytree(cluster_path, os.path.join(build_path, cluster_name))

        # TODO MGnify-MGnify checks

        # Update time
        log({'tot_time': time() - tot_time})
        # Return log
        return log

    @staticmethod
    def iter_clusters(clusters_paths, batch_size=1000, max_clusters=None):
        """Iterate through linclust files

        Args
        linclust_paths (str/list)   Path or list of paths to linclust files
        batch_size (int)            Maximum number of clusters per batch
        max_clusters (int)          Maximum number of clusters to iterate

        Return
        (generator)                 A generator that yelds batch as
                                    tuple(batch index, cluster names)
        """
        # Initialize current cluster index
        cluster_index = 0
        # Initialize current batch index
        batch_index = 0
        # Initialize current batch of clusters
        batch_clusters = list()
        # Loop through each linclust input path
        for linclust_path in get_paths(clusters_paths):
            # Open current linclust path
            with open(linclust_path, 'r') as linclust_file:
                # Loop through each line in current linclust file path
                for linclust_line in linclust_file:
                    # Check if current line is matches expected format
                    match = re.search(r'^(\S+)\s+(\d+)', linclust_line)
                    # Case line does not match expected format: skip
                    if not match: continue
                    # Otherwise, retrieve cluster name and cluster size
                    cluster_name = str(match.group(1))
                    cluster_size = str(match.group(2))
                    # Add current cluster name to batch
                    batch_clusters.append(cluster_name)
                    # Update current cluster index
                    cluster_index += 1
                    # Define batch index
                    batch_index = (cluster_index - 1) // batch_size
                    # Case cluster index has reached maximum size, exit loop
                    if (max_clusters is not None) and (cluster_index >= max_clusters):
                        break
                    # Check if cluster index has not reached batch size
                    elif (cluster_index % batch_size) != 0:
                        # Go to next iteration
                        continue
                    # Case cluster size has reached batch size
                    else:
                        # Yield batch index and list of clusters
                        yield batch_index, batch_clusters
                        # Reset batch of clusters
                        batch_clusters = list()
                # Case cluster index has reached maximum size, exit loop
                if (max_clusters is not None) and (cluster_index >= max_clusters):
                    break
        # In case we reached this point, check for non returned clusters
        if len(batch_clusters) > 0:
            # Yield last batch
            yield batch_index, batch_clusters

    # Wrapper for mgseed pipeline
    def mgseed(self, cluster_names, clusters_dir, min_jobs=1, max_jobs=100):
        # Define log path
        log_path = os.path.join(clusters_dir, 'mgseed.json')
        # Run mgseed pipeline, retrieve log
        self.seed_pipeline(
            cluster_names=cluster_names,
            clusters_dir=batch_path,
            log_path=log_path,
            verbose=True
        )

    # Wrapper for pfbuild pipeline
    def pfbuild(
        self, cluster_names, clusters_dir, target_len=None, target_ds='uniprot',
        log_filename='pfbuild.json', min_jobs=1, max_jobs=100
    ):
        # Initialize list of target datasets paths
        target_paths = list()
        # choose dataset
        if target_ds == 'uniprot':
            # Set uniprot chunks paths
            target_paths = self.uniprot_paths
        elif target_ds == 'mgnify':
            # Set mgnify chunks paths
            target_paths = self.mgnify_paths
        else:
            # Raise error
            raise ValueError('Chosen dataset does not exist')
        # Run pfbuild pipeline
        self.hmm_pipeline(
            cluster_names=cluster_names,
            clusters_dir=clusters_dir,
            target_path=target_paths,
            z_score=target_len,
            e_value=0.01,
            min_jobs=min_jobs,
            max_jobs=max_jobs,
            log_path=os.path.join(clusters_dir, log_filename)
        )

    # Wrapper for mfbuild
    def mfbuild(
        self, cluster_names, clusters_dir, target_len=None, target_ds='mgnify',
        log_filename='mfbuild.json', min_jobs=1, max_jobs=100
    ):
        # Just call pfbuild with different default parameters
        return self.pfbuild(
            cluster_names=cluster_names,
            clusters_dir=clusters_dir,
            target_len=target_len,
            target_ds=target_ds,
            log_filename='pfbuild.json',
            min_jobs=min_jobs,
            max_jobs=max_jobs
        )

    # Get length of target fasta dataset
    def fasta_len(self, target_paths, min_jobs=1, max_jobs=100):
        # Create new client
        client = self.get_client({'cores': 1, 'memory': '2 GB'})
        # Spawn jobs
        client.cluster.adapt(minimum=min_jobs, maximum=max_jobs)
        # Get fasta datasets
        fasta = Fasta.from_list(target_paths)
        # Get each chunk size
        futures = client.map(len, fasta)
        # Retrieve results
        results = client.gather(futures)
        # Shutdown client
        client.shutdown()
        # Compute and return z score
        return = sum(results)
