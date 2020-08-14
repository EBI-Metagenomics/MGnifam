# Dependencies
from src.hmm.hmm import HMM
import os
import re


class Cluster(object):

    # Constructor
    def __init__(
        self, acc='', id='', name='', desc='', auth='', type='', seed='', path='',
        seq_scores=(25.0, 25.0, 25.0), dom_scores=(25.0, 25.0, 25.0)
        # build_method='', search_method=''
    ):
        # Save attributes
        self.acc = acc  # Accession (primary key)
        self.id = id  # Identifier
        self.name = name
        self.desc = desc
        self.auth = auth
        self.type = type
        self.seed = seed
        self.path = path
        # Store (TC, NC, GA) scores for sequences
        self.seq_scores = seq_scores
        # Store (TC, NC, GA) scores for domains
        self.dom_scores = dom_scores
        # # Define build method
        # self.build_method = build_method
        # # Define search method
        # self.search_method = search_method

    def get_path(self, *args):
        return os.path.join(self.path, *args)

    def to_dict(self):
        return {
            'name': self.name,
            'desc': self.desc,
            'auth': self.auth,
            'type': self.type,
            'seed': self.seed,
            'path': self.path,
            'seq_scores': list(self.seq_scores),
            'dom_scores': list(self.dom_scores)
        }

    # Get line format, return default value, cast retrieved value
    @classmethod
    def is_param(cls, line, param='AC', value=r'(\S+)', default='', cast=str):
        # Check if accession matches expected format
        match = re.search(r'^{:s}\s+{:s}'.format(param, value), line)
        # Case line does not match expected format
        if not match:
            # Return default value
            return default

        # Return parameter value
        return cast(match.group(1))

    @classmethod
    def is_acc(cls, line, default=''):
        return cls.is_param(line=line, param='AC', default=default)

    @classmethod
    def is_id(cls, line, default=''):
        return cls.is_param(line=line, param='ID', default=default)

    @classmethod
    def is_desc(cls, line, default=''):
        return cls.is_param(line=line, param='DE', default=default)

    @classmethod
    def is_auth(cls, line, default=''):
        return cls.is_param(line=line, param='AU', default=default)

    @classmethod
    def is_seed(cls, line, default=''):
        return cls.is_param(line=line, param='SE', default=default)

    @classmethod
    def is_type(cls, line, default=''):
        return cls.is_param(line=line, param='TP', default=default)

    @classmethod
    def is_score(cls, line, param='TC', default=[None, None]):
        # Define floating point regex
        is_float = r'[+-]?[0-9]+\.[0-9]+'
        # Retrieve scores
        scores = cls.is_param(
            line=line, param='TC', default='', cast=str,
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
    def is_tc(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.is_score(line, param='TC', default=default)

    @classmethod
    def is_nc(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.is_score(line, param='NC', default=default)

    @classmethod
    def is_ga(cls, line, default=(None, None)):
        # Retrieve scores
        return cls.is_score(line, param='GA', default=default)

    @classmethod
    def from_desc(cls, path, cluster_name='', cluster_path=''):
        # Initialize cluster parameters dictionary
        params = {
            'acc': '', 'id': '', 'desc': '', 'auth': '', 'type': '',
            'name': cluster_name, 'path': cluster_path,
            'seq_scores': [None, None, None],
            'dom_scores': [None, None, None]
        }
        # Open DESC file
        with open(path, 'r') as file:
            # Loop through each line in file
            for line in file:
                # Case line is cluster accession
                params['acc'] = cls.is_accession(line, params['acc'])
                # Case line is id
                params['id'] = cls.is_id(line, params['acc'])
                # Case line is description
                params['desc'] = cls.is_description(line, params['desc'])
                # Case line is author name
                params['auth'] = cls.is_author(line, params['auth'])
                # Case line is cluster type
                params['type'] = cls.is_type(line, params['type'])
                # Case line is seed sequence
                params['seed'] = cls.is_seed(line, params['seed'])
                # Case line is TC item
                params['seq_scores'][0], params['dom_scores'][0] = cls.is_tc(
                    default=(params['seq_scores'][0], params['dom_scores'][0]),
                    line=line,
                )
                # Case line is NC item
                params['seq_scores'][1], params['dom_scores'][1] = cls.is_nc(
                    default=(params['seq_scores'][1], params['dom_scores'][1]),
                    line=line
                )
                # Case line is GA item
                params['seq_scores'][2], params['dom_scores'][2] = cls.is_ga(
                    default=(params['seq_scores'][2], params['dom_scores'][2]),
                    line=line
                )

    @classmethod
    def from_dir(cls, cluster_path):
        # Retrieve directory name as cluster name
        cluster_name = os.path.basename(cluster_path)
        # Define DESC file path
        desc_path = os.path.join(cluster_path, 'DESC')
        # Parse DESC file, retrieve cluster
        cluster = cls.from_desc(desc_path, cluster_name, cluster_path)
        # Return cluster
        return cluster

    def to_desc(self):
        # Open description file
        with open(self.get_path('DESC'), 'w') as file:
            # Write attributes
            file.write('AC   {:s}\n'.format(self.acc))
            file.write('ID   {:s}\n'.format(self.name))
            file.write('DE   {:s}\n'.format(self.desc))
            file.write('AU   {:s}\n'.format(self.auth))
            # file.write('SE   {:s}\n'.format(self.se))
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
            # Write family
            file.write('TP   {:s}'.format(self.family))
