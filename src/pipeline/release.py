# Dependencies
from src.pipeline.pipeline import Pipeline, Log
from src.pipeline.batch import Batch
from src.pipeline.build import Build
from time import time
import sys
import os
import re


class ReleasePipeline(object):
    """ Make new MGnifam release

    This pipeline puts together three pipelines: first, if required, it creates
    the input dataset such that distributed cluster features can be exploited
    at the best (casically chunking the input datasets).

    Then it runs the BATCH pipeline, which makes SEED alignments out of MGnifam
    dataset, applies a threshold over compositional bias, applies automatic SEED
    trimming, makes HMM models and finally checks them against UniProt to ensure
    there is no itersection.

    Next, it runs the BUILD pipeline, which builds the actual MGnifam entry, by
    retrieving sequences from MGnifam and aligning them (through hmmalign) to
    form the new entry alignment, finally providing annotations and database
    management in order to put MGnifam entries in place.
    """

    # Constructor
    def __init__(
        # Pipeline parameters, required to handle job scheduling
        self, cluster_type, cluster_kwargs,
        # Path to datasets (fixed)
        linclust_path, mgnifam_path, uniprot_path,
        # Wether to reformat input dataset
        make_linclust=False, make_mgnifam=False, make_uniprot=False,
        # Compositional bias threshold settings
        comp_bias_threshold=0.2, comp_bias_inclusive=True,
        # Automatic trimming settings
        auto_trim_threshold=0.4, auto_trim_inclusive=True,
        auto_filter_thresold=0.5, auto_filter_inclsive=True,
        # Post trimming settings
        seed_min_width=1, seed_min_height=1,
        # Search against UniProt settings
        uniprot_e_value=0.01, uniprot_z_score=None,
        # Search against MGnifam settings
        mgnifam_e_value=0.01, mgnifam_z_score=None,
        # Environmental variables
        env=os.environ.copy()
    ):
        # TODO Set inner build pipeline
        # TODO Set upload pipeline
        raise NotImplementedError

    def run(self, cluster_names, clusters_path, batch_size=1, max_clusters=1000, verbose=False):
        """ Make new MGnifam release

        First, loops through cluster names, retrieves sequences associated with
        each cluster, eventually discards them according to a bias threshold,
        aligns them to provide SEED alignments, develops HMM from those
        alignments and discards clusters with UniProt intersections.

        Then, loops through remaining clusters to generate a wider multiple
        sequence alignment, by aligning sequneces matching with cluster's
        HMM models throuh hmmalign.

        Afterwards, checks for HMM models overlap with others already in
        MGnifam or with others among the inserted ones, keeping only the
        biggest clusters in term of sequences.
        """
        # Make abstract
        raise NotImplementedError

    @staticmethod
    def iter_cluster_names(path_list, verbose=False):
        # Define set of cluster names
        cluster_names = set()
        # Loop through each file
        for path in path_list:
            # Safely open file
            try:
                # Open file buffer at current path
                file = open(path, 'r')
                # Loop through every line in file
                for line in file:
                    # Match expected line format
                    match = re.search(r'^(\S+)', line)
                    # Case line does not match expected format
                    if not match:
                        # Go to next line
                        continue
                    # Otherwise, save current line's cluster name
                    cluster_names.add(match.group(1))
                # Close file buffer
                file.close()

            except FileNotFoundError:
                # Verbose log
                if verbose:
                    print('Could not open input file', end=' ')
                    print('at path {:s}'.format(path))
                # Raise exception
                raise

            # Return set of cluster names
            return cluster_names
