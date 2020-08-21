from src.pipeline.pipeline import Pipeline
from src.pipeline.build import parse_env_json
from src.hmm.hmmer import Tblout, Domtblout
from src.hmm.hmmpress import HMMPress
from src.hmm.hmmemit import HMMEmit
from src.hmm.hmmscan import HMMScan
from src.sequences import fasta_iter
from src.database import MGnifam as MGDatabase
from src.database import Pfam as PFDatabase
from src.clusters import MGnifam as MGCluster
# from src.clusters import PFam as PFCluster
from subprocess import CalledProcessError
from tempfile import mkstemp
from glob import glob, iglob
from time import time
import argparse
import shutil
import sys
import os


class Load(Pipeline):
    """ Load MGnifam relese into database

    Starting from a build directory, loads each new MGnifam cluster and
    searches for intersections with either previous MGnifam and Pfam release.
    """

    # Constructor
    def __init__(
        self, mgnifam_db, pfam_db,
        hmmpress_cmd=['hmmpress'],
        hmmemit_cmd=['hmmemit'],
        hmmscan_cmd=['hmmscan'],
        env=os.environ
    ):
        # Store databases
        self.mgnifam_db = mgnifam_db
        self.pfam_db = pfam_db
        # Define scripts instances
        self.hmmpress = HMMPress(cmd=hmmpress_cmd, env=env)
        self.hmmemit = HMMEmit(cmd=hmmemit_cmd, env=env)
        self.hmmscan = HMMScan(cmd=hmmscan_cmd, env=env)
        # Store environmental variables
        self.env = env

    # Run pipeline
    def run(self, release_path, clusters_iter, mgnifam_version, verbose=False):
        """ Load clusters to MGnifam database

        First, iterate through each given cluster path and load it into MGnifam
        dataset. Then, chech MGnifam dataset for intersections among clusters.

        Note: either insertions and checks are done inside the same
        transaction. If some error happens, database is rollback to previous
        state.

        Args
        release_path (str)          Path to release directory
        clusters_iter (list)        List of cluster directories paths
        mgnifam_version (str)       String for mgnifam version
        verbose (bool)              Whether to print out verbose log
        """

        # Try to make new release
        try:

            # Connect to database
            self.mgnifam_db.connect()

            # Load each input cluster in database
            self.load_clusters(
                clusters_iter=clusters_iter,
                mgnifam_version=mgnifam_version,
                verbose=verbose
            )

            # Initialize release directory
            self.init_release(clusters_path=release_path)

            # Retrieve MGnifam clusters
            mgnifam_dict, hmm_len, seed_len, align_len = self.retrieve_clusters(
                clusters_path=release_path,
                db=self.mgnifam_db,
                verbose=verbose
            )

            # Initialize sequence scores dict(accession: (TC, NC, GA))
            sequence_scores = dict()
            # Initialize domain scores dicy(accession: (TC, NC, GA))
            domain_scores = dict()
            # Loop through each cluster in MGnifam
            for cluster_acc, cluster in mgnifam_dict.items():
                # Update sequence scores
                sequence_scores.setdefault(cluster_acc, cluster.seq_scores)
                # Update domain scores
                domain_scores.setdefault(cluster_acc, cluster.dom_scores)

            # Check MGnifam clusters against each other
            mgnifam_dict, discard_dict = self.check_mgnifam(
                # Store clusters directories into release folder
                clusters_path=release_path,
                # Dictionary mapping cluster accession to MGnifam cluster
                mgnifam_dict=mgnifam_dict,
                # Dictionary mapping cluster accession to sequence scores
                sequence_scores=sequence_scores,
                # Dictionary mapping cluster accession to domain scores
                domain_scores=domain_scores,
                # Dictionary mapping cluster accession to SEED alignment length
                seed_len=seed_len,
                # Whether to print out verbose output
                verbose=verbose
            )

            # Loop through all discarded clusters
            for cluster in discard_dict.values():
                # Check that cluster folder exists
                if os.path.isdir(cluster.path):
                    # Remove cluster directory
                    shutil.rmtree(cluster.path)

            # Raise exception
            raise Exception()

            # Retrieve Pfam clusters
            pfam_dict, _, _, _ = self.retrieve_clusters(
                clusters_path=release_path,
                db=self.pfam_db,
                verbose=verbose
            )

            # Loop through each cluster in Pfam
            for cluster_acc, cluster in pfam_dict.items():
                # Update sequence scores
                sequence_scores.setdefault(cluster_acc, cluster.seq_scores)
                # Update domain scores
                domain_scores.setdefault(cluster_acc, cluster.dom_scores)

            # check MGnifam clusters against the ones in Pfam
            self.check_pfam(
                # Store clusters directories into release folder
                clusters_path=release_path,
                # Dictionary mapping cluster accession to MGnifam cluster
                mgnifam_dict=mgnifam_dict,
                # Dictionary mapping cluster accession to Pfam cluster
                pfam_dict=pfam_dict,
                # Dictionary mapping cluster accession to sequence scores
                sequence_scores=sequence_scores,
                # Dictionary mapping cluster accession to domain scores
                domain_scores=domain_scores,
                # Whether to print out verbose output
                verbose=verbose
            )

            # Loop through all discarded clusters
            for cluster in discard_dict.values():
                # Kill entry: move it from active MGnifam to dead ones
                self.mgnifam_db.kill_cluster(accession=cluster.accession)
                # Check that cluster folder exists
                if os.path.isdir(cluster.path):
                    # Remove cluster directory
                    shutil.rmtree(cluster.path)

        # If any exception got catched, rollback
        except Exception:
            # Discard changes
            self.mgnifam_db.rollback()
            # Go on raising exception
            raise

        # Otherwise, commit changes
        self.mgnifam_db.conn.commit()

    def load_clusters(self, clusters_iter, mgnifam_version, verbose=False):
        """ Load clusters into MGnifam database

        Args
        clusters_iter (list)    Iterator though cluster directories
        mgnifam_version (str)   MGnifam version to set in given clusters
        verbose (bool)          Whether to print verbose output
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Loading {:d}'.format(len(clusters_iter)), end=' ')
            print('clusters to MGnifam database...')

        # Define next cluster accession
        next_accession = 0
        # Loop through each cluster path
        for cluster_path in clusters_iter:
            # Try loading cluster into database
            try:
                # Parse cluster from directory
                cluster = MGCluster.from_dir(cluster_path)
                # Set creation and update times to now
                cluster.created = cluster.updated = int(time())
                # Case cluster accession is not set
                if not cluster.accession_ and next_accession:
                    # Set next accession
                    cluster.accession = '{:05d}'.format(next_accession)

                # Load cluster to MGnifam database
                cluster.to_database(self.mgnifam_db)

                # Update next available accession
                next_accession = int(cluster.accession_) + 1

            # Catch any exception
            except Exception:
                # Raise any exception
                raise
            # Case value error has been found
            except ValueError as err:
                # Verbose
                if verbose:
                    # Print discarded cluster
                    print('Cluster {:s} discarded:'.format(cluster.name), end=' ')
                    print(str(err).strip())
                # Go to next iteration
                pass

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('...done in {:.0f} seconds'.format(time_end - time_beg))

    def init_release(self, clusters_path):
        """ Initialize release directory structure

        Args
        clusters_path (str)         Path to release directory
        """
        # Check if release directory is already in place
        if not os.path.isdir(clusters_path):
            # Try making release directory
            os.mkdir(clusters_path)

    def retrieve_clusters(self, clusters_path, db, verbose=False):
        """ Retrieve all MGnifam clusters

        Loop through each accession in MGnifam/Pfam dataset, retrieve it and
        make cluster directory. Fill directory with `HMM.model` file.

        Note that db instance must have `get_accessions`, `get_cluster` and
        `make_hmm_model` methods.

        Args
        clusters_path (str)         Path to release directory
        db (Database)               Database instance to use (MGnifam or Pfam)
        verbose (bool)              Whether to print out verbose log

        Return
        (dict)                      Dictionary mapping cluster accession
                                    numbers to clusters objects
        (dict)                      Dictionary mapping cluster accession
                                    numbers to HMM lengths
        (dict)                      Dictionary mapping cluster accession
                                    numbers to SEED alignments length
        (dict)                      Dictionary mapping cluster accession
                                    numbers to ALIGN alignments length
        """
        # # Set database reference
        # db = self.mgnifam_db

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Retrieving MGnifam clusters...')

        # Initialize clusters dictionary
        clusters_dict = dict()
        # Initialize HMM model and alignemnts lengths dictionaries
        hmm_len, seed_len, align_len = dict(), dict(), dict()
        # Retrieve all sequence accessions from MGnifam database
        clusters_acc = db.get_accessions()
        # Loop through all retrieved accessions
        for cluster_acc in clusters_acc:
            # Define cluster path
            cluster_path = os.path.join(clusters_path, cluster_acc)

            # Try making cluster directory
            try:
                # Case cluster  path does not exist
                if not os.path.isdir(cluster_path):
                    # Make cluster directory
                    os.mkdir(cluster_path)

                # Retrieve cluster results and items lengths
                retrieved = db.get_cluster(
                    acccession=cluster_acc
                )
                # Retrieve cluster
                cluster = retrieved[0]

                # # Make DESC file
                # cluster.to_desc()

                # # Define HMM model path
                # model_path = os.path.join(cluster_path, 'HMM.model')

                # Set cluster path
                cluster.path = cluster_path
                # Make HMM model
                db.make_hmm_model(
                    accession=cluster_acc,
                    path=cluster.model_path
                )

                # Store HMM length
                hmm_len.setdefault(cluster_acc, int(retrieved[1]))
                # Store SEED alignment length (number of sequences)
                seed_len.setdefault(cluster_acc, int(retrieved[2]))
                # Store ALIGN alignment length
                align_len.setdefault(cluster_acc, int(retrieved[3]))
                # Append cluster to clusters dictionary
                clusters_dict.setdefault(cluster_acc, cluster)

            # Catch exceptions
            except Exception:
                # Verbose
                if verbose:
                    # Print cluster discarded
                    print('  Could not retrieve {:s} cluster'.format(cluster_acc))
                # Discard cluster directory
                shutil.rmtree(cluster_path)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('...done in {:0f} seconds'.format(time_end -  time_beg))

        # Return clusters dictionary
        return clusters_dict, hmm_len, seed_len, align_len

    def make_consensus(self, clusters_path, verbose=False):
        """ Generate consensus FASTA file

        First, loop through each HMM in mgnifam and run hmmemit
        Then, append fasta sequences retireved from hmmemit to consensus
        fasta file

        Args
        clusters_path (set)         Path to folder holding clusters directories
        verbose (bool)              Whether to print verbose output
        """
        # Define HMM models iterator
        models_iter = glob(os.path.join(clusters_path, '*', 'HMM.model'))

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making consensus for', edn=' ')
            print('{:d} clusters'.format(len(models_iter)), end=' ')
            print('at {:s}...'.format(clusters_path), end=' ')

        # Define all clusters FASTA consensus file
        cons_path = os.path.join(clusters_path, 'CONS.fa')
        # Open all-clusters consensus FASTA file
        with open(cons_path, 'w') as cons_file:
            # Loop through each HMM model
            for model_path in models_iter:
                # Get cluster path
                cluster_path = os.path.dirname(model_path)
                # Define consensus fasta file path
                fasta_path = os.path.join(cluster_path, 'CONS.fa')

                # Just call HMM emit
                self.hmmemit(
                    model_path=model_path,
                    out_path=fasta_path
                )

                # Define consensus fasta sequences
                fasta_sequences = dict()
                # Open output fasta file
                with open(fasta_path, 'r') as fasta_file:
                    # Iterate through each fasta entry in given file
                    for fasta_entry in fasta_iter(fasta_file):
                        # Split entry into header and residues
                        head, resid = tuple(fasta_entry.split('\n'))
                        # Store fasta sequence
                        fasta_sequences.setdefault(head, resid)

                # Get cluster accession (cluster folder is named after accession)
                cluster_acc = os.path.basename(cluster_path)
                # Open fasta file in write mode
                with open(fasta_path, 'w') as fasta_file:
                    # Iterate through each stored fasta tuple(header, residues)
                    for head, resid in fasta_sequences.items():
                        # Write out fasta line(s)
                        fasta_file.write('>' + cluster_acc + '\n' + resid + '\n')

                # Open fasta file in read mode
                with open(fasta_path, 'r') as fasta_file:
                    # Copy model to library
                    shutil.copyfileobj(fasta_file, cons_file)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    def prepare_scan(self, clusters_path, verbose=False):
        """ Prepare files for hmmscan

        First, make single HMM library by joining together HMM models stored
        in each cluster folder. Then, generate `.h3[pmif]` HMM dataset files.

        Args
        clusters_path (set)         Path to folder holding clusters directories
        verbose (bool)              Whether to print verbose output

        Raise
        (FileNotFoundError)         In case some output file has not been found
        (Exception)                 In case generic exception happened, such as
                                    errors during hmmpress execution
        """
        # Define clusters HMM model iterator
        models_iter = glob(os.path.join(clusters_path, '*', 'HMM.model'))

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Making HMM library and database', end=' ')
            print('out of {:d} HMM models'.format(len(models_iter)), end=' ')
            print('at {:s}...'.format(clusters_path), end=' ')

        # Define HMM library path
        library_path = os.path.join(clusters_path, 'HMM.library')
        # Open output file
        with open(library_path, 'w') as library_file:
            # Loop through each HMM model file
            for model_path in models_iter:
                # Open input file
                with open(model_path, 'r') as model_file:
                    # Copy model to library
                    shutil.copyfileobj(model_file, library_file)

        # Try running hmmpress
        try:
            # Run hmmpress script against library
            self.hmmpress(model_path=library_path, overwrite=True)
            # Initialize generated file paths
            h3p_path = next(iglob(clusters_path + '/*.h3p'))
            h3m_path = next(iglob(clusters_path + '/*.h3m'))
            h3i_path = next(iglob(clusters_path + '/*.h3i'))
            h3f_path = next(iglob(clusters_path + '/*.h3f'))

            # Check if all the files are in place
            if not (h3p_path and h3m_path and h3i_path and h3f_path):
                # Raise file not found exception
                raise FileNotFoundError('\n'.join([
                    'some files have not been generated:',
                    '- h3p: {:s}'.format(os.path.basename(h3p_path)),
                    '- h3m: {:s}'.format(os.path.basename(h3m_path)),
                    '- h3i: {:s}'.format(os.path.basename(h3i_path)),
                    '- h3f: {:s}'.format(os.path.basename(h3f_path))
                ]))

        # Intercept process error
        except CalledProcessError as err:
            # Raise exception
            raise Exception('\n'.join([
                'hmmpress script returned {}:'.format(err.returncode),
                err.stderr.strip()
            ]))

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

    def make_scan(self, clusters_path, verbose=False):
        """ Scan MGnifam

        Scan all MGnifam clusters in root directory, using previously
        generated consensus sequences FASTA file and previously generated
        HMM model/library.

        Args
        clusters_path (str)         Path to clusters directory
        verbose (bool)              Whether to print verbose output

        Return
        (str)                       Path to output file
        (str)                       Path to tblout file
        (str)                       Path to domtblout file

        Raise
        (FileNotFoundError)         In case some file has not been found
        (Exception)                 In case generic exception happened, such as
                                    errors during hmmpress execution
        """
        # Define HMM library path
        library_path = os.path.join(clusters_path, 'HMM.library')
        # Initialize FASTA consensus file path
        cons_path = os.path.join(clusters_path, 'CONS.fa')

        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Scanning clusters at', end=' ')
            print('{:s}...'.format(clusters_path), end=' ')

        # Case cannot find HMM library file
        if not os.path.isfile(library_path):
            # Raise exception
            raise FileNotFoundError(
                'could not find HMM library at {:s}'.format(library_path)
            )

        # Case cannot find HMM library file
        if not os.path.isfile(cons_path):
            # Raise exception
            raise FileNotFoundError(
                'could not find consensus FASTA file at {:s}'.format(cons_path)
            )

        # Try executing hmmscan
        try:
            # Define tblout path
            tblout_path = mkstemp(dir=clusters_path, suffix='.tblout')
            # Define domtblout path
            domtblout_path = mkstemp(dir=clusters_path, suffix='.domtblout')
            # Run hmmscan, make .domtblout and .tblout files
            self.hmmscan(
                # Set input HMM model path
                # NOTE: .h3[mifp] files must be in the same directory
                hmm_path=library_path,
                # Set input FASTA file path
                fasta_path=cons_path,
                # # Define output path
                # out_path=out_path,
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

        # Intercept hmmscan execution error
        except CalledProcessError as err:
            # Raise exception
            raise Exception('hmmscan returned {}:\n{}'.format(
                                err.returncode,
                                err.stderr.strip()
                            ))

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} seconds'.format(time_end - time_beg))

        # Return either tblout and domtblout paths
        return tblout_path, domtblout_path

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
                                (keys) to list of sequence hits (values)
        (dict)                  Dictionary mapping cluster accession numbers
                                (keys) to list of domain hits (values)
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0
            # Show execution start
            print('Parsing sequnece and domain hits...', end=' ')

        # Parse tblout
        tblout = Tblout(tblout_path)
        # Retireve hits
        sequence_hits = tblout.get_hits(scores=sequence_scores)

        # Parse domtblout
        domtblout = Domtblout(domtblout_path)
        # Retrieve hits
        domain_hits = domtblout.get_hits(scores=domain_scores)

        # Remove tblout and domtblout temporary files
        os.unlink(tblout_path), os.unlink(domtblout_path)

        # Verbose
        if verbose:
            # Update timers
            time_end = time()
            # Show execution time
            print('done on {.0f} seconds'.format(time_end - time_beg))

        # Return either sequence and domain hits
        return sequence_hits, domain_hits

    def check_mgnifam(self, clusters_path, mgnifam_dict, sequence_scores, domain_scores, seed_len, verbose=False):
        """ Check MGnifam clusters against each other

        Search MGnifam clusters against each other: first, generates consensus
        sequences fasta file for each cluster, which will be merged into a
        single fasta file whose sequence accession number are cluster accession
        numbers themselves, and a single HMM library by merging all HMM models
        retrieved from MGnifam database.

        Then, parses results: if there is a cluster whose SEED alignment
        is lower than half with respect to SEED alignment of another one, and
        those two have at least one significan intersection, then remove the
        former one.

        Args
        clusters_path (str)     Path to directory where clusters will be stored
        mgnifam_dict (dict)     Dictionary mapping MGnifam cluster accession
                                numbers to MGnifam cluster instances
        sequence_scores (dict)  Dictionary associating cluster accession
                                numbers to `(TC, NC, GA)` sequence scores
        domain_scores (dict)    Dictionary associating cluster accession
                                numbers to `(TC, NC, GA)` domain scores
        seed_len (dict)         Dictionary associating cluster accession
                                numbers to their SEED alignment length
        verbose (bool)          Whether to print out verbose log

        Return
        (dict)                  Dictionary mapping MGnifam cluster accession
                                numbers to kept MGnifam cluster instances
        (dict)                  Dictionary mapping MGnifam cluster accession
                                numbers to discarded MGnifam cluster instances
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0

        # Make consensus out of MGnifam dict
        self.make_consensus(clusters_path=clusters_path, verbose=verbose)

        # Prepare hmmscan
        self.prepare_scan(clusters_path=clusters_path, verbose=verbose)

        # Make MGnifam scan
        tblout_path, domtblout_path = self.make_scan(
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Parse MGnifam scan
        sequence_hits, domain_hits = self.parse_scan(
            # Path to sequence hits file
            tblout_path=tblout_path,
            # Path to domain hits file
            domtblout_path=domtblout_path,
            # Significant sequences scores
            sequence_scores=sequence_scores,
            # Significant domain scores
            domain_scores=domain_scores
        )

        # Define discarded clusters dict(accession: cluster)
        discard_dict = dict()
        # Make a set of cluster accessions, ordered according to SEED size
        clusters_acc = sorted(seed_len.keys(), key=seed_len.get)
        # Loop through each cluster accession (query name) in domains hit
        for cluster_acc in clusters_acc:
            # Retrieve hits associated to current cluster accession
            cluster_hits = domain_hits.get(cluster_acc, [])
            # Case current accession has not in sequence hits
            if not sequence_hits.get(cluster_acc, []):
                # Go to next iteration
                continue

            # Case current cluster has been discarded
            if cluster_acc in set(discard_dict.keys()):
                # Skip execution, go to next iteration
                continue

            # Go through each target name (other cluster accession number)
            for hit in cluster_hits:
                # Define target accession
                target_acc = hit.get('target_name')
                # Case target accession equals current cluster accession
                if target_acc == cluster_acc:
                    # Go to next iteration
                    continue

                # Case target cluster has already been discarded
                if target_acc in set(discard_dict.keys()):
                    # Skip execution, go to next iteration
                    continue

                # Define current cluster SEED length
                cluster_seed_len = seed_len.get(cluster_acc, 0)
                # Define target cluster SEED length
                target_seed_len = seed_len.get(target_acc, 0)
                # Check case target cluster SEED aligment is big enough
                if target_seed_len > (cluster_seed_len / 2):
                    # Keep it, go to next iteration
                    continue

                # Retrieve current target cluster
                target_cluster = mgnifam_dict.pop(target_acc)
                # Move target cluster to discarded cluster dictionary
                discard_dict.setdefault(target_acc, target_cluster)

        # Verbose
        if verbose:
            # Define total number of clusters
            num_clusters = len(mgnifam_dict) + len(discard_dict)
            # Show number of checked clusters
            print('Checking {:d} clusters against MGnifam...'.format(num_clusters), end=' ')

            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} sceonds:'.format(time_end - time_beg))

            # Show number of kept and number of discarded clusters
            print('  {:d} clusters have been kept'.format(len(mgnifam_dict)))
            print('  {:d} clusters have beed discarded'.format(len(discard_dict)))

        # Return either kept and discarded dictionaries cluster
        return mgnifam_dict, discard_dict

    def check_pfam(self, clusters_path, mgnifam_dict, pfam_dict, sequence_scores, domain_scores, verbose=False):
        """ Check MGnifam clusters against Pfam ones

        Similar to MGnifam-MGnifam checks, it discards each cluster with at
        least one intersection with Pfam, independently of SEED alignment size.

        Args
        clusters_path (str)     Path to directory where clusters will be stored
        mgnifam_dict (dict)     Dictionary mapping MGnifam cluster accession
                                numbers to MGnifam cluster instances
        pfam_dict (dict)        Dictionary mapping Pfam cluster accession
                                numbers to Pfam cluster instances
        sequence_scores (dict)  Dictionary associating cluster accession
                                numbers to `(TC, NC, GA)` sequence scores
        domain_scores (dict)    Dictionary associating cluster accession
                                numbers to `(TC, NC, GA)` domain scores
        verbose (bool)          Whether to print out verbose log

        Return
        (dict)                  Dictionary mapping MGnifam cluster accession
                                numbers to kept MGnifam cluster instances
        (dict)                  Dictionary mapping MGnifam cluster accession
                                numbers to discarded MGnifam cluster instances
        """
        # Verbose
        if verbose:
            # Initialize timers
            time_beg, time_end = time(), 0.0

        # Make consensus out of MGnifam dict
        self.make_consensus(clusters_path=clusters_path, verbose=verbose)

        # Prepare hmmscan
        self.prepare_scan(clusters_path=clusters_path, verbose=verbose)

        # Make MGnifam scan
        tblout_path, domtblout_path = self.make_scan(
            clusters_path=clusters_path,
            verbose=verbose
        )

        # Parse MGnifam scan
        sequence_hits, domain_hits = self.parse_scan(
            # Path to sequence hits file
            tblout_path=tblout_path,
            # Path to domain hits file
            domtblout_path=domtblout_path,
            # Significant sequences scores
            sequence_scores=sequence_scores,
            # Significant domain scores
            domain_scores=domain_scores
        )

        # Define discarded clusters dict(accession: cluster)
        discard_dict = dict()
        # Define clusters dictionary by merging mgnifam and pfam ones
        clusters_dict = dict(**mgnifam_dict, **pfam_dict)
        # Loop through each cluster accession (query name) in domains hit
        for cluster_acc, cluster in clusters_dict.items():
            # Retrieve hits associated to current cluster accession
            cluster_hits = domain_hits.get(cluster_acc, [])
            # Case current accession has not in sequence hits
            if not sequence_hits.get(cluster_acc, []):
                # Go to next iteration
                continue

            # Case current cluster has been discarded
            if cluster_acc in set(discard_dict.keys()):
                # Skip execution, go to next iteration
                continue

            # Go through each target name (other cluster accession number)
            for hit in cluster_hits:
                # Define target accession
                target_acc = hit.get('target_name')
                # Case target accession equals current cluster accession
                if target_acc == cluster_acc:
                    # Go to next iteration
                    continue

                # Case target cluster has already been discarded
                if target_acc in set(discard_dict.keys()):
                    # Skip execution, go to next iteration
                    continue

                # Check if current cluster is a Pfam one
                cluster_is_pfam = cluster_acc in set(pfam_dict.keys())
                # Check if target cluster is a Pfam one
                target_is_pfam = target_acc in set(pfam_dict.keys())
                # Case current and target clusters are both in Pfam
                if cluster_is_pfam and target_is_pfam:
                    # TODO Raise warning if two Pfam cluster intersect

                    # Skip execution, go to next iteration
                    continue

                # Case the current cluster hits Pfam
                if not cluster_is_pfam:
                    # Remove cluster from mgnifam dictionary
                    cluster = mgnifam_dict.pop(cluster_acc)

                # Case target cluster hits Pfam
                if not target_is_pfam:
                    # Remove cluster from mgnifam dictionary
                    cluster = mgnifam_dict.pop(target_acc)

                # Move removed MGnifam cluster to discarded clusters dictionary
                discard_dict.setdefault(cluster.accession, cluster)

        # Verbose
        if verbose:
            # Define total number of clusters
            num_clusters = len(mgnifam_dict) + len(discard_dict)
            # Show number of checked clusters
            print('Checking {:d} clusters against Pfam...'.format(num_clusters), end=' ')

            # Update timers
            time_end = time()
            # Show execution time
            print('done in {:.0f} sceonds:'.format(time_end - time_beg))

            # Show number of kept and number of discarded clusters
            print('  {:d} clusters have been kept'.format(len(mgnifam_dict)))
            print('  {:d} clusters have beed discarded'.format(len(discard_dict)))

        # Return either kept and discarded dictionaries cluster
        return mgnifam_dict, discard_dict


# Command line execution
if __name__ == '__main__':

    # Argparse
    parser = argparse.ArgumentParser(description='Build MGnifam clusters')
    # Define path to input clusters files
    parser.add_argument(
        '-i', '--in_path', nargs='+', type=str, default=[],
        help='Path to input cluster directories'
    )
    # Define path to release directory
    parser.add_argument(
        '-o', '--out_path', type=str, required=True,
        help='Path to release directory'
    )
    # Define MGnifam version
    parser.add_argument(
        '-V', '--version', type=str, default='',
        help='MGnifam release name'
    )
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
    # Define hmmscan executable
    parser.add_argument(
        '--hmmscan_cmd', nargs='+', type=str, default=['hmmscan'],
        help='HMMER3 scan executable'
    )
    # Whether to print verbose output
    parser.add_argument(
        '-v', '--verbose', type=int, default=1,
        help='Print verbose output'
    )
    # # Define E-value significance threshold for both MGnifam and Pfam
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

    # Define Pfam database options
    group = parser.add_argument_group('Pfam database options')
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

    # Retrieve arguments
    args = parser.parse_args()

    # Load environmental variables
    env = os.environ.copy()
    # Case environmental variable file is set
    if args.env_path:
        # Update environmental variables using given file file
        env = {**env, **parse_env_json(args.env_path)}

    # Make MGnifam database instance
    mgnifam_db = MGDatabase(
        user=args.mgnifam_user,
        password=args.mgnifam_password,
        host=args.mgnifam_host,
        port=args.mgnifam_port
    )

    # Make Pfam database instance
    pfam_db = PFDatabase(
        user=args.pfam_user,
        password=args.pfam_password,
        host=args.pfam_host,
        port=args.pfam_port
    )

    # Initialize new load pipeline
    pipeline = Load(
        mgnifam_db=mgnifam_db,
        pfam_db=pfam_db,
        hmmpress_cmd=args.hmmpress_cmd,
        hmmemit_cmd=args.hmmemit_cmd,
        hmmscan_cmd=args.hmmscan_cmd
        # mgnifam_e_value=args.e_value
    )

    # Run pipeline
    pipeline(
        release_path=args.out_path,
        clusters_iter=args.in_path,
        mgnifam_version=args.version,
        verbose=args.verbose
    )
