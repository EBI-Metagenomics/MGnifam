# Dependencies
import mysql.connector as dbms
import json
import os
import re


class MGnifam(object):

    # Define database settings
    def __init__(self, user=None, password=None, database=None, host='localhost', port=3309):
        # Store database settings
        self.user = user
        self.password = password
        self.database = database
        self.host = host
        self.port = port
        # Initialize database connection
        self.conn = None

    # Connect to database
    def connect(self):
        # Set inner connection
        self.conn = dbms.connect(
            user=self.user,
            password=self.password,
            database=self.database,
            host=self.host,
            port=self.port
        )

    # Close connection
    def close(self):
        # Close connection
        self.conn.close()
        # Clear connection
        self.conn = None

    # Retrieve database cursor
    def get_cursor(self):
        return self.conn.cursor()

    # Get next available MGDUF id
    def get_next_mgduf(self):
        """ Retrieve next available MGDUF id

        Loop through each entry in active and dead families (concatenated):
        find biggest number and increase it by 1.

        Note that MGDUF ids are formatted as `^MGDUF(\d+)$`

        Return
        (str)       Next available MGDUF id
        """
        # Define "alive" MGnifam entries cursor
        alive_cursor = self.get_cursor()
        # Search for MGnifam ids
        alive_cursor.execute('SELECT mgnifam_id FROM mgnifam')

        # Define "dead" MGnifam entries cursor
        dead_cursor = self.get_cursor()
        # Search for dead MGnifam ids
        dead_cursor.execute('SELECT mgnifam_id FROM dead_mgnifam')

        # Initialize maximum MGDUF id
        mgduf_max = -1
        # Define list of all results cursors
        cursor_list = [alive_cursor, dead_cursor]
        # Loop through each query result
        for cursor in cursor_list:
            # Loop through each retrieved row
            for row in cursor:
                # Retrieve current MGDUF
                mgduf_str = row[0]
                # Check if current MGDUF matches expected format
                match = re.search(r'^MGDUF(\d+)$', mgduf_str)
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
        return 'MGDUF{:d}'.format(mgduf_max + 1)


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
    mgnifam = MGnifam(**auth_dict)
    # Connect to MySQL database
    mgnifam.connect()

    # Retrieve database cursor
    cursor = mgnifam.get_cursor()
    # Execute query
    cursor.execute('SELECT mgnifam_id FROM mgnifam')
    # Loop through each retrieved result
    for i, row in enumerate(cursor):
        # Print row
        print('{:d}-th row:'.format(i+1))
        print(row)

    # Close connection
    mgnifam.close()
