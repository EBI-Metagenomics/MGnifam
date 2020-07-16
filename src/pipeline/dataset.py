# Dependencies
import sys
import os

# Custom dependencies
from src.pipeline.pipeline import Pipeline, Log
from src.dataset import LinClust, Fasta


class DatasetPipeline(Pipeline):
    """Initialize dataset
    Starting from full, high size datasets, splits them in chunks.
    Actually considered dataset are MGnify, UniProt and LinClust results.
    """

    # def __init__(
    #     # Clusters attributes
    #     self, cluster_type, cluster_kwargs,
    #     # Dataset paths
    #     mgnify_path, uniprot_path, linclust_path
    # ):
    #     # Call parent constructor
    #     super().__init__(cluster_type, cluster_kwargs)
    #     # Store MGnify dataset
    #     self.mgnify_path = mgnify_path
    #     self.uniprot_path = uniprot_path
    #     self.linclust_path = linclust_path

    # Execute pipeline
    def run(
        self, log_path='', verbose=False,
        # LinClust dataset parameters
        linclust_path='', linclust_chunk_size=1e07, linclust_out_path='{:04d}',
        # UniProt dataset parameters
        uniprot_path='', uniprot_chunk_size=1e06, uniprot_out_path='{:04d}',
        # MGnify dataset parameters
        mgnify_path='', mgnify_chunk_size=1e06, mgnify_out_path='{:04d}'
    ):

        # Initialize log
        log = Log(log_path=log_path)

        # Instantiate new client
        client = self.get_client()
        # Request at maximum 3 jobs (one for every dataset)
        client.cluster.adapt(minimum=3, maximum=3)

        # Initialize futures dict(dataset_name: future)
        futures = dict()

        # Case linclust path is defined
        if linclust_path:
            # Define original LinClust dataset
            ds_linclust = LinClust(path=linclust_path)
            # Split in chunks
            futures['linclust'] = client.submit(
                # Chunking method
                ds_linclust.to_chunks,
                # Chunks output path,
                chunk_path=linclust_out_path,
                # Size (number of lines) of every cluster
                chunk_size=linclust_chunk_size
            )
            # Add entry to log
            log({
                'linclust': {
                    'input_path': linclust_path,
                    'output_dir': os.path.dirname(linclust_out_path),
                    'chunk_size': linclust_chunk_size,
                    'num_chunks': 0
                }
            })

        # Case uniprot path is defined
        if uniprot_path:
            # Define original UniProt dataset
            ds_uniprot = Fasta(path=uniprot_path)
            # Split in chunks
            futures['uniprot'] = client.submit(
                # Chunking method
                ds_uniprot.to_chunks,
                # Chunks output path,
                chunk_path=uniprot_out_path,
                # Size (number of fasta entries) of every cluster
                chunk_size=uniprot_chunk_size
            )
            # Add entry to log
            log({
                'uniprot': {
                    'input_path': uniprot_path,
                    'output_dir': os.path.dirname(uniprot_out_path),
                    'chunk_size': uniprot_chunk_size,
                    'num_chunks': 0
                }
            })

        # Case mgnify path is defined
        if mgnify_path:
            # Define original UniProt dataset
            ds_mgnify = Fasta(path=mgnify_path)
            # Split in chunks
            futures['mgnify'] = client.submit(
                # Chunking method
                ds_mgnify.to_chunks,
                # Chunks output path,
                chunk_path=mgnify_out_path,
                # Size (number of fasta entries) of every cluster
                chunk_size=mgnify_chunk_size
            )
            # Add entry to log
            log({
                'mgnify': {
                    'input_path': mgnify_path,
                    'output_dir': os.path.dirname(mgnify_out_path),
                    'chunk_size': mgnify_chunk_size,
                    'num_chunks': 0
                }
            })

        # Wait for every future to finish
        results = client.gather([future for future in futures.values()], 'raise')
        # Add number of chunks to each chunked dataset
        for i in range(len(results)):
            # Get key for i-th future
            key = [*futures.keys()][i]
            # Update number of chunks
            log({
                key: {**log[key], **{
                    'num_chunks': results[i]
                }}
            })

        # # Close cluster connection
        # client.cluster.close()
        # # Remove client
        # del client
