import os


class Cluster(object):

    # Constructor
    def __init__(
        self, name='', desc='', auth='', type='Family', path='',
        seq_scores=(25.0, 25.0, 25.0), dom_scores=(25.0, 25.0, 25.0)
        # build_method='', search_method=''
    ):
        # Save attributes
        self.name = name
        self.desc = desc
        self.auth = auth
        self.type = type
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
            'path': self.path,
            'seq_scores': list(self.seq_scores),
            'dom_scores': list(self.dom_scores)
        }

    def to_desc(self):
        # Open description file
        with open(self.get_path('DESC'), 'w') as file:
            # Write attributes
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
