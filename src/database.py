# Dependencies
from src.hmm.hmm import HMM
from src.msa.msa import MSA
from time import time
import mysql.connector as dbms
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
        self.cursor_ = self.conn_.cursor()

    def close(self):
        # Close connection
        self.conn_.close()
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

    def query(self, *args, **kwargs):
        # Retrieve new cursor
        cursor = self.cursor
        # Use cursor to execute query
        cursor.execute(*args, **kwargs)
        # Return cursor
        return cursor

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

    # Get clusters accessions
    def get_accessions(self):
        """ Retrieve all accession numbers

        Return
        (set)               Set of cluster accession numbers
        """
        # Get cursor
        cursor = self.query("""
        SELECT DISTINCT mgnifam_acc
        FROM mgnifam
        """)

        # Initialize setf cluster accession numbers
        clusters_acc = set()
        # Loop through every result in query
        for row in cursor:
            # Update clusters accession number set
            clusters_acc.add(row[0])

        # Return clusters accession numbers
        return clusters_acc

    # Get next available MGYF accesion
    def get_next_accession(self):
        """ Retrieve next available MGYF accession

        Return
        (str)       Next available MGYF accession
        """
        # Make "alive" MGnifam entries cursor
        alive_cursor = self.query('SELECT mgnifam_acc FROM mgnifam')
        # Make "dead" MGnifam entries cursor
        dead_cursor = self.query('SELECT mgnifam_acc FROM dead_mgnifam')
        # Initialize maximum MGDUF id
        mgyf_num, mgyf_max = 0, 0
        # Define list of all results cursors
        cursor_list = [alive_cursor, dead_cursor]
        # Loop through each query result
        for cursor in cursor_list:
            # Loop through each retrieved row
            for row in cursor:
                # Retrieve current MGDUF
                mgyf_str = row[0]
                # Check if current MGDUF matches expected format
                match = re.search(r'^MGYF(\d+)', mgyf_str)
                # Case retrieved string does not match expected format
                if not match:
                    continue  # Skip to next row
                # Retrieve MGDUF number from string
                mgyf_num = int(match.group(1))
                # Case current MGDUF number is greater than current maximum
                if mgyf_num > mgyf_max:
                    # Update current maximum MGYF number
                    mgyf_max = mgyf_num

        # Return maximum MGDUF number + 1
        return 'MGYF{:05d}'.format(mgyf_max + 1)

    # Get next available MGDUF id
    def get_next_id(self):
        """ Retrieve next available MGDUF id

        Return
        (str)       Next available MGDUF id
        """
        # Make "alive" MGnifam entries cursor
        alive_cursor = self.query('SELECT mgnifam_id FROM mgnifam')
        # Make "dead" MGnifam entries cursor
        dead_cursor = self.query('SELECT mgnifam_id FROM dead_mgnifam')
        # Initialize maximum MGDUF id
        mgduf_num, mgduf_max = 0, 0
        # Define list of all results cursors
        cursor_list = [alive_cursor, dead_cursor]
        # Loop through each query result
        for cursor in cursor_list:
            # Loop through each retrieved row
            for row in cursor:
                # Retrieve current MGDUF
                mgduf_str = row[0]
                # Check if current MGDUF matches expected format
                match = re.search(r'^MGDUF(\d+)', mgduf_str)
                # Case retrieved string does not match expected format
                if not match:
                    continue  # Skip to next row
                # Retrieve MGDUF number from string
                mgduf_num = int(match.group(1))
                # Case current MGDUF number is greater than current maximum
                if mgduf_num > mgduf_max:
                    # Update current maximum MGDUF number
                    mgduf_max = mgduf_num

        # Return maximum MGDUF number + 1
        return 'MGDUF{:04d}'.format(mgduf_max + 1)

    # Make MGnifam HMM library
    def make_hmmlib(self, accessions, path):
        """ Retrieve HMM entries in MGnifam as single HMM model

        Args
        accessions (iterable)   List of clusters accessions used to select
                                clusters whose HMM must be retrieved
        path (str)              Path where output HMM model must be stored

        Return
        (HMM)                   HMM model instance representing the whole
                                retrieved library
        """
        # Retrieve MGnifam database cursor
        cursor = self.cursor
        # Parse accessions to string
        accessions = ', '.join([str(acc) for acc in accessions])
        # Execute join
        cursor.execute("""
            SELECT mgnifam_acc, hmm
            FROM mgnifam_HMM
            WHERE mgnifam_acc IN ({0:s})
        """.format(accessions))

        # Initialize single HMM
        hmmlib = HMM(path=path, length=0)
        # Open new output file
        with open(path, 'w') as file:
            # For each row, retrieve HMM model
            for row in cursor:
                # Make new HMM
                hmm = HMM.from_string(row[1])
                # Update HMM library length
                hmmlib.length = max(hmmlib.length, hmm.length)
                # Write HMM string to file
                file.write(row[1])
        # Retrun HMM library
        return hmmlib

    # Make MGnifam entry HMM
    def make_hmm(self, accession, path):
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
        # Retrieve MGnifam database cursor
        cursor = self.cursor
        # Execute join
        cursor.execute("""
            SELECT mgnifam_acc, hmm
            FROM mgnifam_HMM
            WHERE mgnifam_acc = ({0:s})
        """.format(accession))

        # Check results number
        if not cursor:
            raise KeyError(' '.join([
                'Error: unable to retrieve HMM',
                'for accession {:s}'.fromat(accession)
            ]))

        # Open new output file
        with open(path, 'w') as file:
            # Retrieve HMM model
            hmm_string = next(cursor)
            # Store HMM model
            file.write(hmm_string)

        # Define HMM instance
        hmm = HMM.from_string(hmm_string)
        # Update path
        hmm.path = path
        # Return HMM instance
        return hmm

    # Load MGnifam cluster into database
    def load_cluster(self, cluster):
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
        # Case cluster ACC is not set: retrieve next one
        if not cluster.acc:
            # Retrieve next available accession
            cluster.acc = self.get_next_accession()

        # Case cluster ID is not set: retrieve next one
        if not cluster.id:
            # Retrieve next available id
            cluster.id = self.get_next_id()

        # TODO Check if given id is already taken

        # TODO Check if given accession is already taken

        # Check that HMM model file exists
        if not os.path.isfile(cluster.get_path('HMM')):
            # Raise new exception
            raise FileNotFoundError(' '.join([
                'Error: unable to find HMM model',
                'for cluster {:s}'.format(cluster.name),
                'at {:s}'.format(cluster.path)
            ]))
        # Load HMM model
        hmm = HMM.from_file(cluster.get_path('HMM'))

        # check that SEED file exists
        if not os.path.isfile(cluster.get_path('SEED')):
            # Raise new exception
            raise FileNotFoundError(' '.join([
                'Error: unable to find SEED alignment',
                'for cluster {:s}'.format(cluster.name),
                'at {:s}'.format(cluster.path)
            ]))
        # Load SEED alignment
        seed = MSA.from_aln(cluster.get_path('SEED'))

        # Check that ALIGN file exists
        if not os.path.isfile(cluster.get_path('ALIGN')):
            # Raise new exception
            raise FileNotFoundError(' '.join([
                'Error: unable to find ALIGN alignment',
                'for cluster {:s}'.format(cluster.name),
                'at {:s}'.format(cluster.path)
            ]))
        # Load ALIGN alignment
        align = MSA.from_aln(cluster.get_path('ALIGN'))

        # Set standard description
        cluster.desc = 'Protein of unknown function ({:s})'.format(cluster.id)

        # Define update an creation times
        created = updated = time()

        # Make a query to insert cluster values
        self.query("""
        INSERT INTO mgnifam (
            mgnifam_acc, mgnifam_id, description,
            author, seed_source, type,
            sequence_GA, domain_GA,
            sequence_TC, domain_TC,
            sequence_NC, domain_NC,
            model_length, num_seed, num_full,
            updated, created
        )
        VALUES (
            '{mgnifam_acc:s}', '{mgnifam_id:str}', '{mgnifam_de:s}',
            '{mgnifam_au:s}', '{mgnifam_se:s}', '{mgnifam_tp:s}',
            {sequence_ga:.02f}, {domain_ga:.02f},
            {sequence_tc:.02f}, {domain_tc:.02f},
            {sequence_nc:.02f}, {domain_nc:.02f},
            {hmm_len:d}, {seed_len:d}, {align_len:d},
            {updated:d}, {created:d}
        );
        """.format(
            mgnifam_acc=cluster.acc, mgnifam_id=cluster.id,
            mgnifam_de=cluster.desc, mgnifam_tp=cluster.type,
            sequence_ga=cluster.seq_scores[2], domain_ga=cluster.dom_scores[2],
            sequence_tc=cluster.seq_scores[0], domain_tc=cluster.dom_scores[0],
            sequence_nc=cluster.seq_scores[1], domain_nc=cluster.dom_scores[1],
            hmm_len=hmm.length, seed_len=seed.aln.shape[0], align_len=align.aln.shape[0],
            updated=updated, created=created
        ))

        # Return updated cluster
        return cluster

    # Load HMM into database
    def load_hmm(self, accession, path):
        # Initialize model string
        model = ''
        # Open HMM model file
        with open(path, 'r') as file:
            # Read and save HMM model file content
            model = file.read()

        # Make query
        self.query("""
            INSERT INTO mgnifam_hmm(mgnifam_acc, hmm)
            VALUES ('{acc:s}', '{hmm:s}')
        """.format(
            acc=accession, hmm=model
        ))

    # Load SEED alignment
    def load_seed_alignment(self, accession, path):

        # Initialize SEED alignment string
        seed = ''
        # Open SEED alignment file
        with open(path, 'r') as file:
            # Read and save SEED file content
            seed = file.read()

        # Define query string
        query = """
            INSERT INTO mgnifam_seed (mgnifam_acc, seed)
            VALUES ('{acc:s}', '{aln:s}')
        """
        # Format query string
        query = query.format(acc=accession, aln=seed)
        # Make query
        self.query(query)


class Pfam(Database):

    def __init__(self, *args, **kwargs):
        raise NotImplementedError


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
