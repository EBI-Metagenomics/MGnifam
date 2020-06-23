"""Handle sequences reading and writing"""

# Dependencies
from tqdm import tqdm
import tempfile
import os
import re


class Fasta(object):

    @staticmethod
    def read(in_file):
        """Make ieterator over fasta entries

        Args
        in_file (file)  Input fasta file object

        Return
        (iterator)  Allows to iterate over each entry in given fasta file
        """
        # Initialize header and residues (both strings)
        header, residues = '', ''
        # Go through each line in file
        for line in in_file:
            # Check if line matches fasta header
            match_header = re.search(r'^>', line)
            # Case current line is not header
            if not match_header:
                # Check if there is a valid accession number set
                if header:
                    # Store residues in current line
                    residues += re.sub(r'[^a-zA-Z]+', '', line)
            # Case current line is header
            else:
                # Define current header and residues as previous ones
                prev_header, prev_residues = header, residues
                # Set new header and reset residues
                header, residues = re.sub(r'[\n\r]*$', '', line), ''
                # Case previous header was a valid one
                if prev_header and prev_residues:
                    # Return previous header and residues as a single string
                    yield prev_header + '\n' + prev_residues
        # Return last entry, if any
        if header and residues:
            yield header + '\n' + residues


# Unit test
if __name__ == '__main__':

    # Initialize new temporary file
    fasta_file = tempfile.NamedTemporaryFile(delete=False)
    # Fill sample file
    with open(fasta_file.name, 'w') as ff:
        ff.write('\n'.join([
            '>This',
            'AGSGDAGSGASGGDHAGASGASHASJHDHGDHASHASHS',
            '>is',
            'JDHSGSVAHCHSJSVCSSBSJDSJSNSJSJSKNSSJSJD',
            'ABSHAJSBHAJASJASBSAJAJSAKASJ',
            '>a',
            'HDHABAAJSDBSBSCIFUEYRTXBSNCJJGDDHYENSKS',
            'ANSNANASJASJA',
            'AJSASNASASJJCBBXSGYFSFDTFSCBYH',
            '>fasta',
            'MANSJAKSNDJEUEJCNBCGEGDBSDJSDHDGEVDEHSH'
        ]))

    # Read file
    with open(fasta_file.name, 'r') as ff:
        # Loop through every entry
        for entry in tqdm(Fasta.read(ff)):
            # Print single entry
            print(entry)

    # Delete sample file
    os.remove(fasta_file.name)
