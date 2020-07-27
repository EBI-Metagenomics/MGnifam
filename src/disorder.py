# Dependencies
from tqdm import tqdm
import numpy as np
import subprocess
import tempfile
import time
import sys
import os
import re

class MobiDbLite(object):

    # Constructor
    def __init__(self, cmd=['mobidb_lite.py'], env=os.environ.copy()):
        # Save attributes
        self.cmd = cmd  # Path to MobiDB Lite executable
        self.env = env  # Environmental variables

    # Run method
    def run(self, sequences, verbose=True):
        """Run MobiDB Lite disorder predictor
        Takes as input a list of sequences and returns a lis of boolean (1/0)
        values, where an item has value 1 if the residue in the same position is
        predicted to be disordered, 0 otherwise.

        Args
        sequences (dict)        Dictionary whose keys are sequence accession
                                number and values are fasta entries
        verbose (bool)          Wether to show verbose log or not

        Return
        (dict(str: list))       List of boolean (1/0) lists
        """
        # Define temporary input file
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        # Define temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Initialize residues dict(acc: residues)
        residues = dict()
        # Fill residues dictionary
        for acc in sequences.keys():
            # Split sequence in header and residues
            h, r = tuple(sequences[acc].split('\n'))
            # Store residues
            residues[acc] = list(r)

        # Verbose log
        if verbose:
            print('Generating input fasta file...')
        # Make input fasta file
        with open(fasta_file.name, 'w') as ff:
            # Loop through each input sequence
            for acc, seq in tqdm(sequences.items(), disable=(not verbose)):
                # Define new header
                header = '>' + str(acc)
                # Write out fasta entry
                ff.write(header + '\n' + ''.join(residues[acc]) + '\n')

        # Verbose log
        if verbose:
            'Executing MobiDB Lite predictor...'
        # Initialize timers
        time_beg, time_end = time.time(), None
        # Run MobiDB Lite command
        cmd_out = subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # stdout/stderr encoding
            env=self.env,  # Environmental variables
            # Executable arguments
            args=[
                *self.cmd,
                '-o', out_file.name,  # Path to output file
                fasta_file.name  # Path to input fasta file
            ]
        )
        # Update timers
        time_end = time.time()
        # Verbose log
        if verbose:
            print('Command line output:', cmd_out) # Debug
            print('Took {:.0f} seconds to predict disordered regions over {:d} sequences'.format(
                time_end - time_beg,  # Time took to run
                len(sequences)  # Number of input sequences
            ))

        # Verbose log
        if verbose:
            print('Fetching MobiDB Lite results...')
        # Initialize empty dict(acc: residues positions)
        disorder = dict()
        # Read output file
        with open(out_file.name, 'r') as of:
            # Read through each line of MobiDB Lite
            for line in tqdm(of, disable=(not verbose)):
                # Check if line is matching expected format
                match = re.search(r'^(\S+)\s+(\d+)\s+(\d+)', line)
                # Case line does not matches format
                if not match:
                    continue  # Skip iteration
                # Retrieve sequence accession
                acc = match.group(1)
                # Retrieve disordered region boundaries
                dis_beg, dis_end = int(match.group(2)), int(match.group(3))
                # Add current sequence accession entry if any
                disorder.setdefault(acc, [])
                # Add disordered region as list of ones and zeroes
                disorder.get(acc).append([
                    # Set 1 if resiude is disordered, 0 otherwise
                    int(i in set(range(dis_beg-1, dis_end)))
                    # Loop through each position in residues list
                    for i in range(len(residues[acc]))
                ])

        # Add non-predicted entries (as arrays of zeroes only)
        for acc in set(sequences.keys()) - set(disorder.keys()):
            # Get #residues for current accession
            n = len(residues[acc])
            # Create an array of zeroes with length #resdiues
            disorder[acc] = [[0] * n]

        # Delete temporary files
        os.remove(fasta_file.name)
        os.remove(out_file.name)

        # Return disordered regions
        return disorder


def compute_comp_bias(sequences):
    """Computes computational bias
    Compositional bias is just the rate of disordered residues over total
    number of residues in a set of sequences.

    Args
    sequences (dict(str: list))     Dictionary associating accession numbers
                                    to list of predictions

    Return
    (float)                         Compositional bias (rate of disordered
                                    residues over total number of residues)
    """
    # Initialize total number of residues
    num_residues = 0
    # Initialize number of disordered residues
    num_disorder = 0

    # Go through each sequence accession number
    for acc in sequences.keys():
        # Define residues matrix (#regions, #residues)
        seq = np.array(sequences[acc], dtype=np.int)
        # Define number of regions and number of residues
        n, m = seq.shape
        # Update number of disordered residues
        num_disorder += seq.sum()
        # Update number of total residues
        num_residues += m

    # Compute compositional bias
    comp_bias = num_disorder / num_residues
    # Return compositional bias
    return comp_bias


def compute_threshold(sequences, threshold=0, inclusive=False):
    """Apply threshold over disorder predictions
    Defines wether a region i actually disordered or not by applying a
    threshold over number of positive predictions (1) for each residue.

    Args
    sequences (dict(str: list))     Dictionary associating accession numbers
                                    to list of predictions
    threshold (int/str)             Number of potitive predictions needed to
                                    consider a residue actually positive
    inclusive (bool)                Wether to use an inclusive threshold
                                    (greater or equal) or an exclusive
                                    threshold (just greater)

    Return
    (dict(str: list))               Dictionary associating accession numbers
                                    to thresholded residues lists (not a
                                    list of lists anymore)
    """
    # Initialize an empty output sequences dictionary
    disorder = dict()
    # Loop through each input sequences accession numbers
    for acc in sequences.keys():
        # Retrieve current accession predictions (summed along rows)
        pred = np.array(sequences[acc], dtype=np.int).sum(axis=0)

        # Case all predictions must be positive to return positive
        if threshold == 'all':
            # Define number of rows and columns
            n, m = pred.shape
            # Apply threshold over sum over columns
            pred = (pred >= n)

        # Case threshold is inclusive
        elif inclusive:
            # Apply threshold to sum over columns
            pred = (pred >= threshold)

        # Case threshold is exclusive
        else:
            # Apply threshold to sum over columns
            pred = (pred > threshold)

        # Turn boolean vector into integer 0/1 vector
        pred = pred.astype(np.int)

        # # Case no residue has been positively predicted
        # if not pred.sum():
        #     continue  # Skip iteration

        # Otherwise, save thresholded prediction as list
        disorder[acc] = [pred.tolist()]

    # Return thresholded disordered regions
    return disorder


# Unit test
if __name__ == '__main__':

    # Get environmental variables
    env = {**os.environ.copy(), **{
        'PATH': ':'.join([
            "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/make",
            "/nfs/production/xfam/pfam/software/bin",
            "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/mgnifam",
            "/usr/lib64/qt-3.3/bin",
            "/ebi/lsf/ebi/ppm/10.2/bin",
            "/ebi/lsf/ebi/ppm/10.2/linux2.6-glibc2.3-x86_64/bin",
            "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/etc",
            "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/bin",
            "/usr/lpp/mmfs/bin",
            "/usr/local/bin",
            "/usr/bin",
            "/usr/local/sbin",
            "/usr/sbin",
            "/bin",
            "/usr/bin",
            "/homes/dclementel/bin"
        ]),
        'PYTHONPATH': ':'.join([
            '/ebi/sp/pro1/interpro/python-modules/lib64/python',
            '/ebi/sp/pro1/interpro/python-modules/lib/python'
        ])
    }}

    # We know that cluster MGYP001051398202 has been identified as biased
    # cluster, so we use its members to test MobiDB Lite python wrapper

    # Define list of test fasta entries
    sequences = {
        "MGYP000851679260": ">MGYP000851679260 PL=10 UP=0 BIOMES=1000000000000 LEN=23 CR=0\nMYEIKKLLTLILIAALTMTPLPV",
        "MGYP001032022277": ">MGYP001032022277 PL=10 UP=0 BIOMES=0000000000001 LEN=109 CR=0\nGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
        "MGYP001036154158": ">MGYP001036154158 PL=10 UP=0 BIOMES=0000000000001 LEN=224 CR=0\nLVACGGGSGAGDNNTPSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDYITFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPKNENNSHPGAWDVEFPDSNNPAYLRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
        "MGYP001038205783": ">MGYP001038205783 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=1\nMKKLLALILAAALALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
        "MGYP000111223802": ">MGYP000111223802 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=1\nMKKTLALILTAALTLSLAACGGSNGTENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001043698415": ">MGYP001043698415 PL=10 UP=0 BIOMES=0000000000001 LEN=115 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHKTLSFIINHA",
        "MGYP001044920212": ">MGYP001044920212 PL=11 UP=0 BIOMES=0000000000001 LEN=232 CR=0\nKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYI",
        "MGYP001045957942": ">MGYP001045957942 PL=01 UP=0 BIOMES=0000000000001 LEN=126 CR=0\nVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001051398202": ">MGYP001051398202 PL=00 UP=0 BIOMES=0000000000001 LEN=240 CR=1\nMKKLLALMLAAALALSLVACGGGGAEDNNTTPSTGNGDTTGTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCITFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPASDYTDCIMEQAYLVTDRYEYTGIPKNENNSHPGAWDVEFPDSDNPEYLRLVYFDDSIDVSQYVGKEITFSAKVMDGDYFESAFVDYIDAIIIE",
        "MGYP001051968473": ">MGYP001051968473 PL=01 UP=0 BIOMES=0000000000001 LEN=52 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTASGGVEDSTPSMTK",
        "MGYP001058115450": ">MGYP001058115450 PL=11 UP=0 BIOMES=0000000000001 LEN=176 CR=0\nDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSQDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNKNDSFPGAWNVEFPNAANRQYLRLVYFDDT",
        "MGYP001061453680": ">MGYP001061453680 PL=01 UP=0 BIOMES=0000000000001 LEN=219 CR=0\nACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
        "MGYP001064539059": ">MGYP001064539059 PL=10 UP=0 BIOMES=0000000000001 LEN=202 CR=0\nMKKLLALMLAAALALSLVACGGGNGAGDTNTSSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDD",
        "MGYP001065068340": ">MGYP001065068340 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALVLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001067122102": ">MGYP001067122102 PL=01 UP=0 BIOMES=0000000000001 LEN=219 CR=0\nLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSIDVGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001067638354": ">MGYP001067638354 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001069194978": ">MGYP001069194978 PL=10 UP=0 BIOMES=0000000000001 LEN=224 CR=0\nALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
        "MGYP001071076496": ">MGYP001071076496 PL=10 UP=0 BIOMES=0000000000001 LEN=55 CR=1\nMKKLLALILAAALALSLVACGGGGGVGDNNTPSTGNGDTTGTDTPSGGEDSTQTA",
        "MGYP001075906476": ">MGYP001075906476 PL=00 UP=0 BIOMES=0000000000001 LEN=240 CR=1\nMKKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
        "MGYP001075906285": ">MGYP001075906285 PL=11 UP=0 BIOMES=0000000000001 LEN=176 CR=1\nSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNYRRLVYFDDSIDVSQYVGKEITFSAKA",
        "MGYP001081128765": ">MGYP001081128765 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001084165602": ">MGYP001084165602 PL=10 UP=0 BIOMES=0000000000001 LEN=62 CR=1\nMKKLLALILAAALALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSIPSGQYLNYSRIESGQ",
        "MGYP001086225311": ">MGYP001086225311 PL=10 UP=0 BIOMES=0000000000001 LEN=74 CR=0\nMKKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILS",
        "MGYP000355976475": ">MGYP000355976475 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001108197041": ">MGYP001108197041 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=0\nMKKLLALILAAALALSLVACGGGSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVILYLPSQDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNKNDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
        "MGYP001113077337": ">MGYP001113077337 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001114448909": ">MGYP001114448909 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALTLSLAACGGSNGTENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLIADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
        "MGYP001125457394": ">MGYP001125457394 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=0\nMKKLLALMLAAALALSLVACGGGNGAGDTNTSSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
        "MGYP001130951151": ">MGYP001130951151 PL=01 UP=0 BIOMES=0000000000001 LEN=141 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITT",
        "MGYP001143804030": ">MGYP001143804030 PL=01 UP=0 BIOMES=0000000000001 LEN=48 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTP"
    }
    # Debug
    print('There are {:d} input sequences'.format(len(sequences)))

    # Make new instance of MobiDB Lite
    mdb_lite = MobiDbLite(env=env)
    # Run MobiDB Lite predictor, get disordered regions
    disorder = mdb_lite.run(sequences, verbose=True)
    # Debug
    print('Predictions:')
    # print(json.dumps(disorder, indent=1))
    for acc in disorder.keys():
        print('  Sequence accession number:', acc)
        for i in range(len(disorder[acc])):
            print('    [' + ', '.join([str(r) for r in disorder[acc][i]]) + ']')

    # Compute threshold version of predicted disorder
    disorder = compute_threshold(disorder, threshold=1, inclusive=1)
    # Debug
    print('Thresholded predictions:')
    for acc in disorder.keys():
        print('  Sequence accession number:', acc)
        for i in range(len(disorder[acc])):
            print('    [' + ', '.join([str(r) for r in disorder[acc][i]]) + ']')

    # Compute compositional bias
    comp_bias = compute_comp_bias(disorder)
    # Debug
    print('Compositional bias: {:.02f}'.format(comp_bias))
