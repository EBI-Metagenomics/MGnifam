"""Handle sequences reading and writing"""

# Dependencies
from tqdm import tqdm
from tempfile import NamedTemporaryFile
import sys
import os
import re


def fasta_iter(iterable):
    """Generate fasta sequences iterator

    Args
    iterable (object)   Iterable object returning a fasta file row when a
                        query is made

    Return
    (generator)         Iterates through each fasta sequence, returning the
                        sequence header and sequence residues separated by a
                        newline character
    """
    # Initialize header and residues (both strings)
    header, residues = '', ''
    # Go through each line in file
    for line in iterable:
        # Check if line matches fasta header
        match_header = re.search(r'^>', line)
        # Case current line is not header
        if not match_header:
            # Check if there is a valid accession number set
            if header:
                # Store residues in current line
                residues += re.sub(r'[^a-zA-Z\-]+', '', line)
        # Case current line is header
        else:
            # Define current header and residues as previous ones
            prev_header, prev_residues = header, residues
            # Set new header and reset residues
            header, residues = re.sub(r'[\n\r]+', '', line), ''
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
    fasta_file = NamedTemporaryFile(delete=False)
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
