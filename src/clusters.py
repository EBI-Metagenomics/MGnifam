# Dependencies
from src.hmm.hmmer import Domtblout
from src.hmm.hmm import HMM
from src.msa.msa import MSA
import os
import re


class Cluster(object):

    # Constructor
    def __init__(
        self, accession='', id='', description='', author='', source='', type='', path='',
        seq_scores=(25.0, 25.0, 25.0), dom_scores=(25.0, 25.0, 25.0),
        created=None, updated=None
    ):
        # Store path to file
        self.path = path
        # Save attributes
        self.accession = accession  # Accession (primary key)
        self.id = id  # Identifier
        self.description = description  # Cluster description
        self.author = author  # Name of the cluster annotator
        self.type = type  # Type of the cluster
        self.created = created  # When cluster was created
        self.updated = updated  # When clusters was last updated
        # Store (TC, NC, GA) scores for sequences
        self.seq_scores = seq_scores
        # Store (TC, NC, GA) scores for domains
        self.dom_scores = dom_scores

    @property
    def desc_path(self):
        return os.path.join(self.path, 'DESC')

    @property
    def hmm_path(self):
        return os.path.join(self.path, 'HMM.model')

    @property
    def seed_path(self):
        return os.path.join(self.path, 'SEED.aln')

    @property
    def align_path(self):
        return os.path.join(self.path, 'ALIGN.aln')

    @property
    def hits_path(self):
        return os.path.join(self.path, 'HITS.tsv')


class MGnifam(Cluster):

    # Constructor
    def __init__(
        self, accession='', id='', name='', version='', author='', path='',
        seq_scores=(25.0, 25.0, 25.0), dom_scores=(25.0, 25.0, 25.0),
        created=None, updated=None
    ):
        # Store path to file
        self.path = path
        # Save attributes
        self.accession = accession  # Accession (primary key)
        self.id = id  # Identifier
        self.name = name  # Name of the cluster (SEED sequence)
        self.author = author  # Cluster annotator name
        self.version = version  # MGnifam version
        self.created = created  # When cluster was created
        self.updated = updated  # When clusters was last updated
        # Store (TC, NC, GA) scores for sequences
        self.seq_scores = seq_scores
        # Store (TC, NC, GA) scores for domains
        self.dom_scores = dom_scores

    @property
    def accession(self):
        return 'MGYF{:s}'.format(self.accession_)

    @accession.setter
    def accession(self, accession):
        self.accession_ = re.sub(r'^MGYF', '', accession)

    @property
    def id(self):
        return 'MGDUF{:s}'.format(self.id_)

    @id.setter
    def id(self, id):
        self.id_ = re.sub(r'^MGDUF', '', id)

    @property
    def description(self):
        return 'Protein of unknown function ({:s})'.format(self.id)

    @property
    def type(self):
        return 'Family'

    @property
    def source(self):
        return '{0:s} ({1:s})'.format(self.name, self.version)

    @classmethod
    def from_dir(cls, cluster_path):
        # Retrieve directory name as cluster name
        cluster_name = os.path.basename(cluster_path)
        # Define DESC file path
        desc_path = os.path.join(cluster_path, 'DESC')
        # Parse DESC file, retrieve cluster and return it
        return cls.from_desc(desc_path, cluster_name, cluster_path)

    def to_desc(self):
        # Open description file
        with open(self.desc_path, 'w') as file:
            # Write attributes
            file.write('AC   {:s}\n'.format(self.accession))
            file.write('ID   {:s}\n'.format(self.id))
            file.write('DE   {:s}\n'.format(self.description))
            file.write('AU   {:s}\n'.format(self.author))
            # Write gathering threshold (GA) for sequence, domain pair
            file.write('GA   {seq_ga:.02f} {dom_ga:.02f};\n'.format(
                seq_ga=self.seq_scores[2],
                dom_ga=self.dom_scores[2]
            ))
            # Write upper threshold (TC) for sequence, domain pair
            file.write('TC   {seq_tc:.02f} {dom_tc:.02f};\n'.format(
                seq_tc=self.seq_scores[0],
                dom_tc=self.dom_scores[0]
            ))
            # Write lower threshold (NC) for sequence, domain pair
            file.write('NC   {seq_nc:.02f} {dom_nc:.02f};\n'.format(
                seq_nc=self.seq_scores[1],
                dom_nc=self.dom_scores[1]
            ))
            # Wirte source
            file.write('SE   {:s}'.format(self.source))
            # Write family
            file.write('TP   {:s}'.format(self.type))

    # Parse to dictionary
    def to_dict(self):
        return {
            'accession': self.accession,
            'id': self.id,
            'description': self.description,
            'author': self.author,
            'type': self.type,
            'source': self.source,
            'seq_scores': list(self.seq_scores),
            'dom_scores': list(self.dom_scores),
            'name': self.name,
            'path': self.path,
            'created': self.created,
            'updated': self.updated
        }

    @classmethod
    def parse_param(cls, line, param='AC', value=r'(\S+)', default='', cast=str):
        # Check if accession matches expected format
        match = re.search(r'^{:s}\s+{:s}'.format(param, value), line)
        # Case line does not match expected format
        if not match:
            # Return default value
            return default

        # Return parameter value
        return cast(match.group(1))

    @classmethod
    def parse_accession(cls, line, default=''):
        return cls.parse_param(line=line, param='AC', default=default)

    @classmethod
    def parse_id(cls, line, default=''):
        return cls.parse_param(line=line, param='ID', default=default)

    @classmethod
    def is_desc(cls, line, default=''):
        return cls.parse_param(line=line, param='DE', default=default)

    @classmethod
    def parse_author(cls, line, default=''):
        return cls.parse_param(line=line, param='AU', default=default)

    @classmethod
    def parse_source(cls, line, default=''):
        return cls.parse_param(line=line, param='SE', default=default)

    @classmethod
    def parse_name(cls, line, default=''):
        # Return full parameter string `<name> (<version>)`
        param = cls.parse_source(line=line, default='')

        # Case line is empty
        if param:
            # Try matching name
            match = re.search(r'^(\S+)', param)

            # Case name is not matched
            if match:
                # Otherwise, return retrieved name
                return str(match.group(1))

        # Return default value
        return default

    @classmethod
    def parse_version(cls, line, default=''):
        # Return full parameter string `<name> (<version>)`
        param = cls.parse_source(line=line, default='')

        # Case line is not empty
        if param:
            # Try matching version
            match = re.search(r'\((.*)\)$', param)

            # Case version is matched
            if match:
                # Otherwise, return retrieved name
                return str(match.group(1))

        # Return default value
        return default

    @classmethod
    def parse_score(cls, line, param='TC', default=[None, None]):
        # Define floating point regex
        is_float = r'[+-]?[0-9]+\.[0-9]+'
        # Retrieve scores
        scores = cls.parse_param(
            line=line, param=param, default='', cast=str,
            value='({0:s}\s+{0:s})[;]?'.format(is_float)
        )
        # Case no score has been found
        if not scores:
            # Return default ones
            return default

        # Split scores in sequence and domain score
        seq_score, dom_score = tuple([float(s) for s in re.split(r'\s+', scores)])
        # Return scores
        return seq_score, dom_score

    @classmethod
    def parse_tc(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.parse_score(line, param='TC', default=default)

    @classmethod
    def parse_nc(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.parse_score(line, param='NC', default=default)

    @classmethod
    def parse_ga(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.parse_score(line, param='GA', default=default)

    @classmethod
    def from_desc(cls, path, cluster_name='', cluster_path=''):
        # Initialize cluster parameters dictionary
        params = {
            'accession': '', 'id': '', 'author': '', 'version': '',
            'name': cluster_name, 'path': cluster_path,
            'seq_scores': [None, None, None],
            'dom_scores': [None, None, None]
        }

        # Open DESC file
        with open(path, 'r') as file:
            # Loop through each line in file
            for line in file:
                # Case line is cluster accession
                params['accession'] = cls.parse_accession(line, params['accession'])
                # Case line is id
                params['id'] = cls.parse_id(line, params['id'])
                # Case line is author name
                params['author'] = cls.parse_author(line, params['author'])
                # Retrieve name and version out pf source line
                params['name'] = cls.parse_name(line, params['name'])
                params['version'] = cls.parse_version(line, params['version'])
                # Case line is TC item
                params['seq_scores'][0], params['dom_scores'][0] = cls.parse_tc(
                    default=(params['seq_scores'][0], params['dom_scores'][0]),
                    line=line,
                )
                # Case line is NC item
                params['seq_scores'][1], params['dom_scores'][1] = cls.parse_nc(
                    default=(params['seq_scores'][1], params['dom_scores'][1]),
                    line=line
                )
                # Case line is GA item
                params['seq_scores'][2], params['dom_scores'][2] = cls.parse_ga(
                    default=(params['seq_scores'][2], params['dom_scores'][2]),
                    line=line
                )
        # Return cluster instance
        return cls(**params)

    # Load cluster to MGnifam database
    def to_database(self, db):
        """ Load cluster to MGnifam database

        Args
        db (database.MGnifam)               MGnifam database instance, must be
                                            authenticated

        Raise
        (ValueError)                        In some fields are not correct
        (FileNotFoundError)                 In case one of the required files
                                            has not been found
        (OSError)                           In case one of the parsed files is
                                            not correctly formatted
        """
        # Retrieve next available accession number
        next_accession = db.get_next_accession()
        # Case accession is not set
        if not self.accession_:
            # Set accession number
            self.accession = next_accession

        # Check accession (next accession must be greater than current)
        if int(self.accession_) < int(re.sub('^MGYF', '', next_accession)):
            # Raise exception
            raise ValueError(' '.join([
                'Could not load cluster into database:',
                'current accession is {:s}, while'.format(self.accession),
                'next available accession {:s}'.format(next_accession)
            ]))

        # Retrieve next available id
        next_id = db.get_next_id()
        # Case id is not set
        if not self.id_:
            # Set accession number
            self.id = next_id

        # Check id (next id must be greater than current one)
        if int(self.id_) < int(re.sub('^MGDUF', '', next_id)):
            # Raise exception
            raise ValueError(' '.join([
                'Could not load cluster into database:',
                'current id is {:s}, while'.format(self.id),
                'next available id {:s}'.format(next_id)
            ]))

        # Load HMM model file
        hmm_model = HMM.from_file(self.hmm_path)
        # Check HMM model name
        if not hmm_model.name:
            # Raise exception
            raise ValueError('HMM model has no name set')
        # Check HMM model length
        if not hmm_model.length:
            # Raise exception
            raise ValueError('HMM model length is not valid')

        # Load SEED alignment
        seed_msa = MSA.from_aln(self.seed_path)
        # Check SEED alignment file
        if seed_msa.is_empty:
            # Raise exception
            raise ValueError('SEED alignment shape is not valid')

        # Check ALIGN alignment file
        align_msa = MSA.from_aln(self.align_path)
        # Check ALIGN file
        if align_msa.is_empty:
            # Raise exception
            raise ValueError('ALIGN alignment shape is not valid')

        # Check HITS.tsv file
        hits = Domtblout.hits_from_tsv(self.hits_path)
        # # Case domain hits list is empty
        # if not len(hits):
        #     # Raise exception
        #     raise ValueError('HITS list shape is empty')

        # Load cluster into database
        db.load_cluster(
            cluster=self,
            hmm_model=hmm_model,
            seed_msa=seed_msa,
            align_msa=align_msa
        )

        # Load HMM model into database
        db.load_hmm_model(self.accession, path=self.hmm_path)

        # Load SEED alignment into database
        db.load_seed_msa(self.accession, path=self.seed_path)

        # Load SEED hits into database
        db.load_seed_hits(self.accession, path=self.seed_path)

        # Load ALIGN hits into database
        db.load_align_hits(self.accession, path=self.hits_path)

        # Return updated cluster
        return self


class Pfam(Cluster):

    @property
    def accession(self):
        return 'PF{:s}'.format(self.accession_)

    @accession.setter
    def accession(self, accession):
        self.accession_ = re.sub(r'^PF', '', accession)
