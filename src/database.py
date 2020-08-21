# Dependencies
from src.clusters import MGnifam as MGCluster
from src.clusters import Pfam as PFCluster
from src.hmm.hmmer import Domtblout
from src.hmm.hmm import HMM
from src.msa.msa import MSA
import mysql.connector as dbms
from datetime import datetime
from time import time
import json
import os
import re


class Database(object):

    # Constructor
    def __init__(self, user='', password='', database='', host='', port=3309, autocommit=False):
        # Store database settings
        self.user = user
        self.password = password
        self.database = database
        self.host = host
        self.port = port
        # Initialize autocommit
        self.autocommit_ = bool(autocommit)
        # Initialize database connection
        self.conn_ = None
        # Initialize database cursor
        self.cursor_ = None

    def connect(self):
        # Set inner connection
        self.conn_ = dbms.connect(
            user=self.user,
            password=self.password,
            database=self.database,
            host=self.host,
            port=self.port
        )
        # Define autocommit
        self.autocommit = False
        # Define cursor
        self.cursor_ = self.conn.cursor()

    def close(self):
        # Close connection
        self.conn.close()
        # Clear cursor
        self.cursor_ = None
        # Clear connection
        self.conn_ = None

    def __enter__(self):
        # Connect database
        self.connect()
        # Retrun database instance
        return self

    def __close__(self, type, value, traceback):
        # Close connection
        self.close()
        # Eventually raise exception
        return type

    def execute(self, *args, **kwargs):
        # Use cursor to execute query
        self.cursor.execute(*args, **kwargs)
        # Return cursor
        return self.cursor

    @property
    def conn(self):
        return self.conn_

    @property
    def cursor(self):
        return self.cursor_

    @property
    def autocommit(self):
        return self.conn.autocommit

    @autocommit.setter
    def autocommit(self, autocommit):
        # Update autocommit flag
        self.autocommit_ = bool(autocommit)
        # Update connection
        self.conn.autocommit = bool(autocommit)

    def rollback(self):
        self.conn.rollback()


class MGnifam(Database):

    # Constructor
    def __init__(self, user='', password='', host='', port=3309, autocommit=False):
        # Parsent super constructor
        super().__init__(user=user, password=password, host=host, port=port, database='mgnifam', autocommit=autocommit)

    # Get clusters accessions
    def get_accessions(self):
        """ Retrieve all accession numbers from MGnifam

        Return
        (set)               Set of cluster accession numbers
        """
        # Define query
        query = 'SELECT DISTINCT mgnifam_acc FROM mgnifam'
        # Initialize set of cluster accession numbers
        accessions = [row[0] for row in self.execute(query)]
        # Return clusters accession numbers
        return set(accessions)

    # Get next available MGYF accesion
    def get_next_accession(self):
        """ Retrieve next available MGYF accession

        Return
        (str)       Next available MGYF accession
        """
        # Execute query
        cursor = self.execute(""" SELECT mgnifam_acc FROM mgnifam
                                  UNION
                                  SELECT mgnifam_acc FROM dead_mgnifam; """)
        # Initialize current and maximum MGYF accession number
        mgyf_int, mgyf_max = 0, 0
        # Loop through each retrieved row
        for row in cursor:
            # Retrieve current MGDUF
            mgyf_string = row[0]
            # Check if current MGDUF matches expected format
            match = re.search(r'^MGYF(\d+)', mgyf_string)
            # Case retrieved string matches expected format
            if match:
                # Retrieve MGDUF number from string
                mgyf_int = int(match.group(1))
                # Case current MGDUF number is greater than current maximum
                if mgyf_int > mgyf_max:
                    # Update current maximum MGYF number
                    mgyf_max = mgyf_int

        # Return maximum MGDUF number + 1
        return 'MGYF{:05d}'.format(mgyf_max + 1)

    # Get next available MGDUF id
    def get_next_id(self):
        """ Retrieve next available MGDUF id

        Return
        (str)       Next available MGDUF id
        """
        # Make "alive" MGnifam entries cursor
        cursor = self.execute(""" SELECT mgnifam_id FROM mgnifam
                                  UNION
                                  SELECT mgnifam_id FROM dead_mgnifam; """)
        # Initialize current and maximum MGDUF id
        mgduf_int, mgduf_max = 0, 0
        # Loop through each retrieved row
        for row in cursor:
            # Retrieve current MGDUF
            mgduf_string = row[0]
            # Check if current MGDUF matches expected format
            match = re.search(r'^MGDUF(\d+)', mgduf_string)
            # Case retrieved string matches expected format
            if match:
                # Retrieve MGDUF number from string
                mgduf_int = int(match.group(1))
                # Case current MGDUF number is greater than current maximum
                if mgduf_int > mgduf_max:
                    # Update current maximum MGDUF number
                    mgduf_max = mgduf_int

        # Return maximum MGDUF number + 1
        return 'MGDUF{:04d}'.format(mgduf_max + 1)

    # Make MGnifam entry HMM
    def make_hmm_model(self, accession, path):
        """ Retrieve HMM for a single accession

        Args
        accession (str)         Cluster accession number for which HMM model
                                must be retrieved
        path (str)              Path where HMM model will be stored

        Raise
        (KeyError)              In case there was no HMM associated with given
                                cluster accession number
        (OSError)               In case it was not possible to save HMM model
        """
        # Define query
        query = """ SELECT mgnifam_acc, hmm
                    FROM mgnifam_HMM
                    WHERE mgnifam_acc = '{0:s}'; """
        # Execute query
        cursor = self.execute(query.format(accession))

        # Retrieve first row
        row = next(cursor)
        # Check results number
        if not row:
            # Raise exception
            raise KeyError(' '.join([
                'cannot retrieve HMM',
                'for accession {:s}'.format(accession)
            ]))

        # Open new output file
        with open(path, 'w') as file:
            # Retrieve HMM model
            hmm_string = row[1]
            # Store HMM model
            file.write(hmm_string)

        # Define HMM instance
        hmm = HMM.from_string(hmm_string)
        # Update path
        hmm.path = path
        # Return HMM instance
        return hmm

    # Load HMM model into cluster
    def load_hmm_model(self, accession, path):
        # Open HMM model file
        with open(path, 'r') as file:
            # Read HMM model as string
            model_string = file.read()

            # Define query
            query = """ INSERT INTO mgnifam_HMM (mgnifam_acc, hmm)
                        VALUES ('{:s}', '{:s}') """
            # Set values into query
            query = query.format(accession, model_string)
            # Execute query
            self.execute(query)

    # Load SEED alignment
    def load_seed_msa(self, accession, path):
        """ Load MSA tinto database

        MySQL [mgnifam]> describe mgnifam_seed;
        +-------------+------------+------+-----+---------+-------+
        | Field       | Type       | Null | Key | Default | Extra |
        +-------------+------------+------+-----+---------+-------+
        | mgnifam_acc | varchar(9) | NO   | MUL | NULL    |       |
        | seed        | mediumblob | YES  |     | NULL    |       |
        +-------------+------------+------+-----+---------+-------+
        """
        # Open alignment file
        with open(path, 'r') as file:
            # Read entire file
            msa_string = file.read()

            # Define query
            query = """ INSERT INTO mgnifam_seed (mgnifam_acc, seed)
                        VALUES ('{:s}', '{:s}'); """
            # Set values into query
            query = query.format(accession, msa_string)
            # Execute query
            self.execute(query)

    # Load SEED hits table
    def load_seed_hits(self, accession, path):
        """ Load SEED hits from .tsv file

        Args
        accession (str)     Cluster accession number
        path (str)          Path to SEED alignmnet

        MySQL [mgnifam]> describe mgnifam_reg_seed;
        +----------------+--------------+------+-----+---------+-------+
        | Field          | Type         | Null | Key | Default | Extra |
        +----------------+--------------+------+-----+---------+-------+
        | mgnifam_acc    | varchar(9)   | NO   | MUL | NULL    |       |
        | mgnifamseq_acc | varchar(16)  | NO   | MUL | NULL    |       |
        | seq_start      | mediumint(8) | NO   |     | 0       |       |
        | seq_end        | mediumint(8) | NO   |     | NULL    |       |
        +----------------+--------------+------+-----+---------+-------+
        """
        # Read MSA alignment
        msa = MSA.from_aln(path=path)

        # Initialize array of values
        values = list()
        # Loop through each SEED alignment line
        for i in range(msa.shape[0]):
            # Define current hit values
            value = "('{0:s}', '{1:s}', {2:d}, {3:d})"
            # Set values
            value = value.format(accession, msa.acc[i], msa.beg[i], msa.end[i])
            # Store values for current alignment
            values.append(value)

        # Initialize query
        query = """ INSERT INTO mgnifam_reg_seed (mgnifam_acc, mgnifamseq_acc, seq_start, seq_end)
                    VALUES {:s}; """
        # Add values
        query = query.format(', '.join(values))

        # Case at least one value has been set
        if i > 0:
            # Execute query
            self.execute(query)

    # Load ALIGN hits table
    def load_align_hits(self, accession, path):
        """ Load ALIGN hits from .tsv table

        Args
        accession (str)         Cluster accession number
        path (str)              Path to tabular file containing domain hits

        MySQL [mgnifam]> describe mgnifam_reg_full;
        +-----------------------+-----------------------+------+-----+---------+----------------+
        | Field                 | Type                  | Null | Key | Default | Extra          |
        +-----------------------+-----------------------+------+-----+---------+----------------+
        | auto_mgnifam_reg_full | int(15) unsigned      | NO   | PRI | NULL    | auto_increment |
        | mgnifam_acc           | varchar(9)            | NO   | MUL | NULL    |                |
        | mgnifamseq_acc        | varchar(16)           | NO   | MUL | NULL    |                |
        | seq_start             | mediumint(8)          | NO   |     | 0       |                |
        | seq_end               | mediumint(8)          | NO   |     | 0       |                |
        | ali_start             | mediumint(8) unsigned | NO   |     | NULL    |                |
        | ali_end               | mediumint(8) unsigned | NO   |     | NULL    |                |
        | model_start           | mediumint(8)          | NO   |     | 0       |                |
        | model_end             | mediumint(8)          | NO   |     | 0       |                |
        | domain_bits_score     | double(8,2)           | NO   |     | 0.00    |                |
        | domain_evalue_score   | varchar(15)           | NO   |     | NULL    |                |
        | sequence_bits_score   | double(8,2)           | NO   |     | 0.00    |                |
        | sequence_evalue_score | varchar(15)           | NO   |     | NULL    |                |
        +-----------------------+-----------------------+------+-----+---------+----------------+
        """
        # Read domain hits file from path
        hits = Domtblout.hits_from_tsv(path, sep='\t')

        # Initialize array of values
        values = list()
        # Loop through each hit
        for i, hit in enumerate(hits):
            # Define current hit values
            value = """(
                '{cluster_acc:s}', '{seq_acc:s}', {seq_beg:d}, {seq_end:d},
                {ali_beg:d}, {ali_end:d}, {model_beg:d}, {model_end:d},
                {dom_bits:f}, {dob_eval:f}, {seq_bits:f}, {seq_eval:f}
            )"""

            # Retrieve sequence accession number
            seq_acc = hit['target_name']
            # Retrieve sequence boundaries
            seq_beg, seq_end = hit['envelope_beg'], hit['envelope_end']
            # Retrieve alignment boundaries
            ali_beg, ali_end = hit['alignment_beg'], hit['alignment_end']
            # Retrieve model boundaries
            model_beg, model_end = hit['model_beg'], hit['model_end']
            # Retrieve domain bit-socore and e-value
            dom_bits, dom_eval = hit['bit_score'], hit['e_value']
            # Set sequence bit-score and e-value
            seq_bits = hit['sequence_bit_score']
            seq_eval = hit['sequnece_e_value']
            # Retrieve values for current hit
            values.append(value.format(
                # Set cluster accession number
                cluster_acc=accession,
                # Set sequence accession number and boundaries
                seq_acc=seq_acc, seq_beg=int(seq_beg), seq_end=int(seq_end),
                # Set alignment boundaries
                ali_beg=int(ali_beg), ali_end=int(ali_end),
                # Set model boundaries
                model_beg=int(model_beg), model_end=int(model_end),
                # Set domain bit-score and e-value
                dom_bits=float(dom_bits), dom_eval=float(dom_eval),
                # Set sequence bit-score and e-value
                seq_bits=float(seq_bits), seq_eval=float(seq_eval)
            ))

        # Initialize query
        query = """ INSERT INTO mgnifam_reg_seed (
            mgnifam_acc, mgnifamseq_acc,
            seq_start, seq_end, ali_start, ali_end, model_start, model_end
            domain_bits_score, domain_evalue_score,
            sequence_bits_score, sequence_evalue_score
        )
        VALUES {:s}; """
        # Add values
        query = query.format(', '.join(values))

    # Retrieve cluster
    def get_cluster(self, accession):
        """ Retrieve cluster from database

        Args
        accession (str)         Accession number of the cluster to retrieve

        Return
        (Cluster)               Cluster instance
        (int)                   HMM model length
        (int)                   SEED alignment length (number of sequences)
        (int)                   ALIGN alignment length (number of sequences)

        Raise
        (KeyError)              In case cluster associated to given accession
                                number has not been found
        """
        # Initialize query
        query = """ SELECT mgnifam_acc, mgnifam_id, author,
                           sequence_TC, sequence_NC, sequence_GA,
                           domain_TC, domain_NC, domain_GA,
                           model_length, num_seed, num_full
                    FROM mgnifam
                    WHERE mgnifam_acc = '{:s}'; """
        # Set cluster accession number
        query = query.format(accession)
        # Execute query
        cursor = self.execute(query)

        # Get first retrieved row
        row = next(cursor)
        # Case no row has been retrieved
        if not row:
            # Raise exception
            raise KeyError('couold not find cluster for accession {:s}'.format(accession))

        # Make an instance out of retrieved cluster
        cluster = MGCluster(
            accession=row[0], id=row[1], author=row[2],
            # Set domain scores
            dom_scores=(row[3], row[4], row[5]),
            # Set sequence scores
            seq_scores=(row[6], row[7], row[8])
        )

        # Retrieve HMM model and alignments lengths
        model_len, seed_len, align_len = row[9], row[10], row[11]

        # Return either cluster and lengths
        return cluster, model_len, seed_len, align_len

    # Kill cluster (move to dead MGnifam)
    def kill_cluster(self, accession, author):
        """ Kill MGnifam cluster

        Args
        accession (str)     Accession number of the cluster to kill
        author (str)        Name of the user which killed the cluster

        Raise
        (KeyError)          In case given cluster accession number does
                            not exist
        """
        # Initialize query: retrieve accession and id from MGnifam
        query = """ SELECT mgnifam_acc, mgnifam_id
                    FROM mgnifam
                    WHERE mgnifam_acc = '{:s}'; """
        # Set values into query
        query = query.format(accession)
        # Execute query, retrieve cursor
        cursor = self.execute(query)

        # Get first row
        row = next(cursor)
        # Case first row is empty
        if not row:
            # Raise error
            raise KeyError('cound not find accession {:s}'.format(accession))
        # Otherwise, retrieve either accession and id
        accession, id = tuple(row)

        # Initialize query: insert current accession into dead MGnifam
        query = """ INSERT INTO dead_mgnifam (mgnifam_acc, mgnifam_id, user, killed)
                    VALUES ('{:s}', '{:s}', '{:s}', {:d});"""
        # Set values into query
        query = query.format(accession, id, author, time())
        # Execute query
        self.execute(query)

    # Load MGnifam cluster into database
    def load_cluster(self, cluster, hmm_model, seed_msa, align_msa):
        """ Load cluster into MGnifam database

        First, generates cluster accession and cluster id if those have not been
        previously set, by retrieving them sa previously available ones.
        Then, reads HMM model, SEED and ALIGN alignmentsand retrieves some
        of their characteristics.
        Finally, adds a new entry in `mgnifam` table.

        Args
        cluster (mgnifam.Cluster)       Cluster containing all attributes

        Return
        (mgnifam.Cluster)               Same input cluster, eventually updated
        """
        # Make query string
        query = """
        INSERT INTO mgnifam (
            mgnifam_acc, mgnifam_id, description, author, seed_source, type,
            sequence_TC, sequence_NC, sequence_GA,
            domain_TC, domain_NC, domain_GA,
            model_length, num_seed, num_full,
            buildMethod, searchMethod,
            msv_lambda, msv_mu,
            viterbi_lambda, viterbi_mu,
            forward_lambda, forward_tau,
            updated, created
        )
        VALUES (
            '{mgnifam_acc:s}', '{mgnifam_id:s}', '{mgnifam_de:s}',
            '{mgnifam_au:s}', '{mgnifam_se:s}', '{mgnifam_tp:s}',
            {sequence_tc:.02f}, {sequence_nc:.02f}, {sequence_ga:.02f},
            {domain_tc:.02f}, {domain_nc:.02f}, {domain_ga:.02f},
            {hmm_len:d}, {seed_len:d}, {align_len:d},
            '{build_method:s}', '{search_method:s}',
            {msv_lambda:.02f}, {msv_mu:.02f},
            {viterbi_lambda:.02f}, {viterbi_mu:.02f},
            {forward_lambda:.02f}, {forward_tau:.02f},
            '{updated:s}', '{created:s}'
        ) """

        # Define current created and updated time
        updated = created = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        # Define sequence scores
        sequence_tc, sequence_nc, sequence_ga = tuple(cluster.seq_scores)
        # Define domain scores
        domain_tc, domain_nc, domain_ga = tuple(cluster.dom_scores)
        # Define HMM model length
        hmm_len = len(hmm_model)
        # Defne SEED alignment length
        seed_len = seed_msa.shape[0]
        # Define ALIGN alignment length
        align_len = align_msa.shape[0]
        # Add values to query string
        query = query.format(
            mgnifam_acc=cluster.accession, mgnifam_id=cluster.id,
            mgnifam_de=cluster.description, mgnifam_au=cluster.author,
            mgnifam_se=cluster.source, mgnifam_tp=cluster.type,
            sequence_tc=float(sequence_tc), domain_tc=float(domain_tc),
            sequence_nc=float(sequence_nc), domain_nc=float(domain_nc),
            sequence_ga=float(sequence_ga), domain_ga=float(domain_ga),
            hmm_len=hmm_len, seed_len=seed_len, align_len=align_len,
            build_method='', search_method='',
            msv_lambda=0.0, msv_mu=0.0,
            viterbi_lambda=0.0, viterbi_mu=0.0,
            forward_lambda=0.0, forward_mu=0.0, forward_tau=0.0,
            updated=updated, created=created
        )
        # Execute query
        self.execute(query)


class Pfam(Database):

    # Constructor
    def __init__(self, user='', password='', host='', port=3309, autocommit=False):
        # Parsent super constructor
        super().__init__(user=user, password=password, host=host, port=port, database='pfam_live', autocommit=autocommit)

    def get_accessions(self):
        """ Retrieve all accession numbers from Pfam

        Return
        (set)               Set of cluster accession numbers

        MySQL [pfam_live]> describe pfamA;
        +----------------------+--------------+------+-----+-------------------+-------+
        | Field                | Type         | Null | Key | Default           | Extra |
        +----------------------+--------------+------+-----+-------------------+-------+
        | pfamA_acc            | varchar(7)   | NO   | PRI | NULL              |       |
        | pfamA_id             | varchar(16)  | NO   | UNI | NULL              |       |
        | previous_id          | tinytext     | YES  |     | NULL              |       |
        | description          | varchar(100) | NO   |     | NULL              |       |
        | deposited_by         | varchar(100) | NO   |     | anon              |       |
        | seed_source          | tinytext     | NO   |     | NULL              |       |
        | type                 | varchar(30)  | NO   | MUL | NULL              |       |
        | comment              | longtext     | YES  |     | NULL              |       |
        | sequence_GA          | double(8,2)  | NO   |     | NULL              |       |
        | domain_GA            | double(8,2)  | NO   |     | NULL              |       |
        | sequence_TC          | double(8,2)  | NO   |     | NULL              |       |
        | domain_TC            | double(8,2)  | NO   |     | NULL              |       |
        | sequence_NC          | double(8,2)  | NO   |     | NULL              |       |
        | domain_NC            | double(8,2)  | NO   |     | NULL              |       |
        | buildMethod          | tinytext     | NO   |     | NULL              |       |
        | model_length         | mediumint(8) | NO   |     | NULL              |       |
        | searchMethod         | tinytext     | NO   |     | NULL              |       |
        | msv_lambda           | double(8,2)  | NO   |     | NULL              |       |
        | msv_mu               | double(8,2)  | NO   |     | NULL              |       |
        | viterbi_lambda       | double(8,2)  | NO   |     | NULL              |       |
        | viterbi_mu           | double(8,2)  | NO   |     | NULL              |       |
        | forward_lambda       | double(8,2)  | NO   |     | NULL              |       |
        | forward_tau          | double(8,2)  | NO   |     | NULL              |       |
        | num_seed             | int(10)      | YES  |     | NULL              |       |
        | num_full             | int(10)      | YES  |     | NULL              |       |
        | updated              | timestamp    | NO   |     | CURRENT_TIMESTAMP |       |
        | created              | datetime     | YES  |     | NULL              |       |
        | version              | smallint(5)  | YES  |     | NULL              |       |
        | number_archs         | int(8)       | YES  |     | NULL              |       |
        | number_species       | int(8)       | YES  |     | NULL              |       |
        | number_structures    | int(8)       | YES  |     | NULL              |       |
        | number_ncbi          | int(8)       | YES  |     | NULL              |       |
        | number_meta          | int(8)       | YES  |     | NULL              |       |
        | average_length       | double(6,2)  | YES  |     | NULL              |       |
        | percentage_id        | int(3)       | YES  |     | NULL              |       |
        | average_coverage     | double(6,2)  | YES  |     | NULL              |       |
        | change_status        | tinytext     | YES  |     | NULL              |       |
        | seed_consensus       | text         | YES  |     | NULL              |       |
        | full_consensus       | text         | YES  |     | NULL              |       |
        | number_shuffled_hits | int(10)      | YES  |     | NULL              |       |
        | number_uniprot       | int(10)      | YES  |     | NULL              |       |
        | rp_seed              | tinyint(1)   | YES  |     | NULL              |       |
        | number_rp15          | int(8)       | YES  |     | NULL              |       |
        | number_rp35          | int(8)       | YES  |     | NULL              |       |
        | number_rp55          | int(8)       | YES  |     | NULL              |       |
        | number_rp75          | int(8)       | YES  |     | NULL              |       |
        +----------------------+--------------+------+-----+-------------------+-------+
        """
        # Define query
        query = 'SELECT DISTINCT pfamA_acc FROM pfamA'
        # Initialize set of cluster accession numbers
        accessions = [row[0] for row in self.execute(query)]
        # Return clusters accession numbers
        return set(accessions)

    def make_hmm_model(self, accession, path):
        """ Retrieve HMM for a single accession

        Args
        accession (str)         Cluster accession number for which HMM model
                                must be retrieved
        path (str)              Path where HMM model will be stored

        Raise
        (KeyError)              In case there was no HMM associated with given
                                cluster accession number
        (OSError)               In case it was not possible to save HMM model

        MySQL [pfam_live]> describe pfamA_HMM;
        +-----------+------------+------+-----+---------+-------+
        | Field     | Type       | Null | Key | Default | Extra |
        +-----------+------------+------+-----+---------+-------+
        | pfamA_acc | varchar(7) | NO   | MUL | NULL    |       |
        | hmm       | mediumblob | YES  |     | NULL    |       |
        | logo      | mediumblob | YES  |     | NULL    |       |
        +-----------+------------+------+-----+---------+-------+
        """
        # Define query
        query = """ SELECT pfamA_acc, hmm
                    FROM pfamA_HMM
                    WHERE pfamA_acc = '{0:s}'; """
        # Execute query
        cursor = self.execute(query.format(accession))

        # Check results number
        if not cursor:
            # Raise exception
            raise KeyError(' '.join([
                'could not retrieve HMM',
                'for accession {:s}'.format(accession)
            ]))

        # Open new output file
        with open(path, 'w') as file:
            # Retrieve HMM model
            hmm_string = next(cursor)[1]
            # Store HMM model
            file.write(hmm_string)

        # Define HMM instance
        hmm = HMM.from_string(hmm_string)
        # Update path
        hmm.path = path
        # Return HMM instance
        return hmm

    # Retrieve cluster
    def get_cluster(self, accession):
        """ Retrieve cluster from database

        Args
        accession (str)         Accession number of the cluster to retrieve

        Return
        (Cluster)               Cluster instance
        (int)                   HMM model length
        (int)                   SEED alignment length (number of sequences)
        (int)                   ALIGN alignment length (number of sequences)

        Raise
        (KeyError)              In case cluster associated to given accession
                                number has not been found
        """
        # Initialize query
        query = """ SELECT pfamA_acc, pfamA_id, description, deposited_by, type,
                           sequence_TC, sequence_NC, sequence_GA,
                           domain_TC, domain_NC, domain_GA,
                           model_length, num_seed, num_full
                    FROM mgnifam
                    WHERE mgnifam_acc = '{:s}'; """
        # Set cluster accession number
        query = query.format(accession)
        # Execute query
        cursor = self.execute(query)

        # Get first retrieved row
        row = next(cursor)
        # Case no row has been retrieved
        if not row:
            # Raise exception
            raise KeyError('couold not find cluster for accession {:s}'.format(accession))

        # Make an instance out of retrieved cluster
        cluster = PFCluster(
            # Set cluster attributes
            accession=row[0], id=row[1], description=row[3], author=row[3], type=row[4],
            # Set domain scores
            dom_scores=(row[5], row[6], row[7]),
            # Set sequence scores
            seq_scores=(row[8], row[9], row[10])
        )

        # Retrieve HMM model and alignments lengths
        model_len, seed_len, align_len = row[11], row[12], row[13]

        # Return either cluster and lengths
        return cluster, model_len, seed_len, align_len


# Unit testing
if __name__ == '__main__':

    """
    MySQL [mgnifam]> describe mgnifam;
    +----------------+--------------+------+-----+-------------------+-------+
    | Field          | Type         | Null | Key | Default           | Extra |
    +----------------+--------------+------+-----+-------------------+-------+
    | mgnifam_acc    | varchar(9)   | NO   | PRI |                   |       |
    | mgnifam_id     | varchar(16)  | NO   | UNI | NULL              |       |
    | previous_id    | tinytext     | YES  |     | NULL              |       |
    | description    | varchar(100) | NO   |     | NULL              |       |
    | author         | tinytext     | YES  |     | NULL              |       |
    | deposited_by   | varchar(100) | NO   |     | anon              |       |
    | seed_source    | tinytext     | NO   |     | NULL              |       |
    | type           | varchar(30)  | NO   |     | NULL              |       |
    | comment       mgnifam_ids | longtext     | YES  |     | NULL              |       |
    | sequence_GA    | double(8,2)  | NO   |     | NULL              |       |
    | domain_GA      | double(8,2)  | NO   |     | NULL              |       |
    | sequence_TC    | double(8,2)  | NO   |     | NULL              |       |
    | domain_TC      | double(8,2)  | NO   |     | NULL              |       |
    | sequence_NC    | double(8,2)  | NO   |     | NULL              |       |
    | domain_NC      | double(8,2)  | NO   |     | NULL              |       |
    | buildMethod    | tinytext     | NO   |     | NULL              |       |
    | model_length   | mediumint(8) | NO   |     | NULL              |       |
    | searchMethod   | tinytext     | NO   |     | NULL              |       |
    | msv_lambda     | double(8,2)  | NO   |     | NULL              |       |
    | msv_mu         | double(8,2)  | NO   |     | NULL              |       |
    | viterbi_lambda | double(8,2)  | NO   |     | NULL              |       |
    | viterbi_mu     | double(8,2)  | NO   |     | NULL              |       |
    | forward_lambda | double(8,2)  | NO   |     | NULL              |       |
    | forward_tau    | double(8,2)  | NO   |     | NULL              |       |
    | num_seed       | int(10)      | YES  |     | NULL              |       |
    | num_full       | int(10)      | YES  |     | NULL              |       |
    | updated        | timestamp    | NO   |     | CURRENT_TIMESTAMP |       |
    | created        | datetime     | NO   |     | CURRENT_TIMESTAMP |       |
    +----------------+--------------+------+-----+-------------------+-------+

    MySQL [mgnifam]> describe dead_mgnifam;
    +-------------+-------------+------+-----+-------------------+-------+
    | Field       | Type        | Null | Key | Default           | Extra |
    +-------------+-------------+------+-----+-------------------+-------+
    | mgnifam_acc | varchar(9)  | YES  | UNI | NULL              |       |
    | mgnifam_id  | varchar(40) | NO   |     | NULL              |       |
    | comment     | mediumtext  | YES  |     | NULL              |       |
    | forward_to  | varchar(7)  | YES  |     | NULL              |       |
    | user        | varchar(10) | NO   |     | anon              |       |
    | killed      | timestamp   | NO   |     | CURRENT_TIMESTAMP |       |
    +-------------+-------------+------+-----+-------------------+-------+

    MySQL [mgnifam]> describe mgnifam_HMM;
    +-------------+------------+------+-----+---------+-------+
    | Field       | Type       | Null | Key | Default | Extra |
    +-------------+------------+------+-----+---------+-------+
    | mgnifam_acc | varchar(9) | NO   | MUL | NULL    |       |
    | hmm         | mediumblob | YES  |     | NULL    |       |
    +-------------+------------+------+-----+---------+-------+
    """

    # Define ROOT path
    ROOT_PATH = os.path.dirname(__file__) + '/..'
    # Define TEMP path
    TEMP_PATH = ROOT_PATH + '/tmp'

    # Initialize database authentication dictionary
    auth_dict = dict()
    # Load database credentials dictionary (at tmp/auth.json)
    with open(TEMP_PATH + '/auth.json', 'r') as auth_file:
        # Read authentication file
        auth_dict = json.load(auth_file)

    # Define new MGnifam database
    mgnifam_db = MGnifam(**auth_dict)
    # Connect to MySQL database
    mgnifam_db.connect()

    # Disable autocommit: start transaction
    mgnifam_db.autocommit = False

    # Make example query
    cursor = mgnifam_db.query("""
    SELECT mgnifam_id
    FROM mgnifam
    """)

    # Loop through each retrieved result
    for i, row in enumerate(cursor):
        # Print row
        print('{:d}-th row:'.format(i+1))
        print(row)

    # Try inserting new test entry into dead MGnifam
    cursor = mgnifam_db.query("""
    INSERT INTO dead_mgnifam (mgnifam_id, mgnifam_acc)
    VALUES ('MGYPTEST', 'MGYPTEST')
    """)

    # Rollback: abort changes
    mgnifam_db.rollback()

    # Look for test entry
    cursor = mgnifam_db.query("""
    SELECT *
    FROM dead_mgnifam
    WHERE mgnifam_id = 'MGYPTEST'
    """)

    # Assert no entry found
    assert not cursor.rowcount, 'Error: test entry found into database'

    # Close connection
    mgnifam_db.close()
