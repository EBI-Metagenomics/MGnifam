# Dependencies
import numpy as np
import subprocess
import os
import re


class MobiDbLite(object):

    # Constructor
    def __init__(self, cmd=['mobidb_lite.py'], env=os.environ.copy()):
        # Save attributes
        self.cmd = cmd  # Path to MobiDB Lite executable
        self.env = env  # Environmental variables

    # Wrapper for run method
    def __call__(self, *args, **kwargs):
        # Just call run method
        return self.run(*args, **kwargs)

    # Run method
    def run(self, fasta_path, out_path=''):
        """ Run MobiDB Lite predictor

        Takes as input a fasta file, return a predictions .tsv file

        Args
        fasta_path (str)        Path to input file
        out_path (str)          Path to output file if any

        Raise
        (CalledProcessError)    In case it could not successfully run MobiDB
                                Lite predictor.
        (FileNotFoundError)     In case input file does not exist
        """
        # Check that input file exists
        if not os.path.isfile(fasta_path):
            # Raise exception
            raise FileNotFoundError(' '.join([
                'could not find input file',
                'at {:s}'.format(fasta_path)
            ]))

        # Initialize command line arguments
        args = self.cmd
        # Add output file, if any
        args += ['-o', out_path] if out_path else []
        # Add fasta file (mandatory)
        args += [fasta_path]

        # Run script, return result
        return subprocess.run(
            capture_output=True,  # Capture console output
            encoding='utf-8',  # stdout/stderr encoding
            check=True,  # Check that script worked properly
            env=self.env,  # Environmental variables
            args=args
        )

    @staticmethod
    def parse(out_path):
        """ Parse MobiDB Lite output file

        Args
        out_path (str)          Path to disorder predictions output path

        Return
        (disorder.Disorder)     Disordered prediction regions instance
        """
        # Just return disorder predictions instance
        return Disorder(path=out_path)


class Disorder:

    # Constructor
    def __init__(self, path):
        # Save file path
        self.path = path
        # Initialize disordered regions predictions
        self.regions = dict()
        # Open file at given path
        with open(self.path, 'r') as file:
            # Loop through each line in file
            for line in file:
                # Match expected format from file string
                match = re.search(r'^(\S+)\s+(\d+)\s+(\d+)', line)
                # Case expected format is not matches
                if not match:
                    # Skip iteration, go to next line
                    continue

                # Retrieve accession and disordered region boundaries
                acc, beg, end = str(match.group(1)), int(match.group(2)), int(match.group(3))
                # Initialize predicted disordered regions list
                self.regions.setdefault(acc, list())
                # Update predicted disordered regions list
                self.regions[acc] += [(beg, end)]

    # Compute compositional bias over some sequences
    def comp_bias(self, fasta_sequences, threshold=1, inclusive=True):
        """ Compute compositional bias

        Compute compositional bias for given dictionary, mapping sequence
        accession numbers to fasta entries. Threshold is applied over
        single residue majority voting prediction to define whether it must be
        considered disordered or not

        Args
        fasta_sequences (dict)      Dictionary mapping sequence accession
                                    numbers to fasta entries, used to retrieve
                                    total sequneces length
        threshold (int)             How many times a single residue must be
                                    predicted disordered in order to actually
                                    consider it disordered.
        inclusive (bool)            Whether the given threshold must be
                                    considered inclusive or exclusive

        Return
        (float)                     Compositional bias computed over all
                                    given sequences
        """
        # Initialize number of residues and disordered ones in given sequences
        num_residues, num_disorder = 0, 0
        # Loop through each accession in given sequences
        for acc in fasta_sequences:
            # Split current fasta entry into header and residues
            head, resid = tuple(fasta_sequences.get(acc).split('\n'))
            # Compute current number of reisudes
            n = len(resid)

            # Initialize zeroes binary sequences prediction
            sequence_bin = np.zeros(n, dtype=np.int)
            # Loop through predicted disordered regions for current accession
            for beg, end in self.regions.get(acc, []):
                # Update current region in binary sequence
                sequence_bin[beg-1:end] += 1

            # Case threshold is inclusive
            if inclusive:
                # Apply threshold over binary sequence
                sequence_bin = sequence_bin >= threshold
            # Case threshold is exclusive
            if not inclusive:
                # Apply threshold over binary sequence
                sequence_bin = sequence_bin > threshold

            # Update number of residues
            num_residues += n
            # Update number of disordered residues
            num_disorder += np.sum(sequence_bin)

        # Compute and return compositional bias
        return num_disorder / num_residues

    # # Run method
    # def run(self, fasta_sequences, out_path=None, verbose=False):
    #     """Run MobiDB Lite disorder predictor
    #     Takes as input a list of sequences and returns a lis of boolean (1/0)
    #     values, where an item has value 1 if the residue in the same position is
    #     predicted to be disordered, 0 otherwise.
    #
    #     Args
    #     fasta_sequences (dict)  Dictionary whose keys are sequence accession
    #                             numbers and values are fasta entries
    #     out_path (str)          Where to save MobiDB output file
    #     verbose (bool)          Wether to show verbose log or not
    #
    #     Return
    #     (dict)                  Dict associating to accession numbers lists of
    #                             boolean (1/0) lists
    #     """
    #     # Define temporary input file
    #     fasta_path = NamedTemporaryFile(delete=False).name
    #
    #     # Check if output file path has been set or must be made
    #     make_out = (out_path is None)
    #     # If output file must be made
    #     if make_out:
    #         # Define temporary output file
    #         out_path = NamedTemporaryFile(delete=False).name
    #
    #     # Verbose log
    #     if verbose:
    #         print('Generating input fasta file...')
    #     # Make input fasta file
    #     with open(fasta_path, 'w') as fasta_file:
    #         # Loop through each input sequence
    #         for sequence_acc, fasta_sequence in fasta_sequences.items():
    #             # Split fasta entry in header and residues
    #             header, residues = tuple(fasta_sequence.split('\n'))
    #             # Define new header using only sequence accession number
    #             header = '>' + str(sequence_acc)
    #             # Write out fasta entry
    #             fasta_file.write(header + '\n' + ''.join(residues) + '\n')
    #
    #     # Verbose log
    #     if verbose:
    #         'Executing MobiDB Lite predictor...'
    #     # Initialize timers
    #     time_beg, time_end = time.time(), None
    #     # Run MobiDB Lite command
    #     ran = subprocess.run(
    #         capture_output=True,  # Capture console output
    #         encoding='utf-8',  # stdout/stderr encoding
    #         check=True,  # Check that script worked properly
    #         env=self.env,  # Environmental variables
    #         args=[*self.cmd, '-o', out_path, fasta_path]
    #     )
    #     # Update timers
    #     time_end = time.time()
    #
    #     # Verbose log
    #     if verbose:
    #         print('Command line output:', ran)
    #         print('Took {:.0f} seconds to predict disordered regions over {:d} sequences'.format(
    #             time_end - time_beg,  # Time took to run
    #             len(fasta_sequences)  # Number of input sequences
    #         ))
    #
    #     # Verbose log
    #     if verbose:
    #         print('Fetching MobiDB Lite results...')
    #     # Initialize empty dict(acc: residues positions)
    #     disorder = dict()
    #     # Read output file
    #     with open(out_path, 'r') as out_file:
    #         # Read through each line of MobiDB Lite
    #         for line in tqdm(out_file, disable=(not verbose)):
    #             # Check if line is matching expected format
    #             match = re.search(r'^(\S+)\s+(\d+)\s+(\d+)', line)
    #             # Case line does not matches format
    #             if not match:
    #                 continue  # Skip iteration
    #             # Retrieve sequence accession
    #             acc = match.group(1)
    #             # Retrieve disordered region boundaries
    #             dis_beg, dis_end = int(match.group(2)), int(match.group(3))
    #             # Add current sequence accession entry if any
    #             disorder.setdefault(acc, [])
    #             # Split associated fasta sequence in header and residues
    #             header, residues = tuple(fasta_sequences.get(acc).split('\n'))
    #             # Add disordered region as list of ones and zeroes
    #             disorder.get(acc).append([
    #                 # Set 1 if resiude is disordered, 0 otherwise
    #                 int(i in set(range(dis_beg-1, dis_end)))
    #                 # Loop through each position in residues list
    #                 for i in range(len(residues))
    #             ])
    #
    #     # Add non-predicted entries (as arrays of zeroes only)
    #     for acc in set(fasta_sequences.keys()) - set(disorder.keys()):
    #         # Get header and residues for current sequence
    #         header, residues = tuple(fasta_sequences.get(acc).split('\n'))
    #         # Create an array of zeroes with length #resdiues
    #         disorder[acc] = [[0] * len(residues)]
    #
    #     # Delete temporary FASTA file
    #     os.remove(fasta_path)
    #
    #     # Eventually, remove output temporary file
    #     if make_out:
    #         os.remove(out_path)
    #
    #     # Return disordered regions
    #     return disorder

    @staticmethod
    def compute_bias(sequences):
        """Computes computational bias
        Compositional bias is just the rate of disordered residues over total
        number of residues in a set of sequences.

        Args
        sequences (dict)    Dictionary associating accession numbers to list of
                            disorder predictions

        Return
        (float)             Compositional bias (rate of disordered residues over
                            total number of residues)
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

    @staticmethod
    def apply_threshold(sequences, threshold=0, inclusive=False):
        """Apply threshold over disorder predictions
        Defines wether a region is actually disordered or not by applying a
        threshold over number of positive predictions (1) for each residue.

        Args
        sequences (dict)                Dictionary associating accession numbers
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

            # Otherwise, save thresholded prediction as list
            disorder[acc] = [pred.tolist()]

        # Return thresholded disordered regions
        return disorder


# # Unit test
# if __name__ == '__main__':
#
#     # Get environmental variables
#     env = {**os.environ.copy(), **{
#         'PATH': ':'.join([
#             "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/make",
#             "/nfs/production/xfam/pfam/software/bin",
#             "/nfs/production/metagenomics/mgnifams/dclementel/Pfam/PfamScripts/mgnifam",
#             "/usr/lib64/qt-3.3/bin",
#             "/ebi/lsf/ebi/ppm/10.2/bin",
#             "/ebi/lsf/ebi/ppm/10.2/linux2.6-glibc2.3-x86_64/bin",
#             "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/etc",
#             "/ebi/lsf/ebi/10.1/linux3.10-glibc2.17-x86_64/bin",
#             "/usr/lpp/mmfs/bin",
#             "/usr/local/bin",
#             "/usr/bin",
#             "/usr/local/sbin",
#             "/usr/sbin",
#             "/bin",
#             "/usr/bin",
#             "/homes/dclementel/bin"
#         ]),
#         'PYTHONPATH': ':'.join([
#             '/ebi/sp/pro1/interpro/python-modules/lib64/python',
#             '/ebi/sp/pro1/interpro/python-modules/lib/python'
#         ])
#     }}
#
#     # We know that cluster MGYP001051398202 has been identified as biased
#     # cluster, so we use its members to test MobiDB Lite python wrapper
#
#     # Define list of test fasta entries
#     fasta_sequences = {
#         "MGYP000851679260": ">MGYP000851679260 PL=10 UP=0 BIOMES=1000000000000 LEN=23 CR=0\nMYEIKKLLTLILIAALTMTPLPV",
#         "MGYP001032022277": ">MGYP001032022277 PL=10 UP=0 BIOMES=0000000000001 LEN=109 CR=0\nGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
#         "MGYP001036154158": ">MGYP001036154158 PL=10 UP=0 BIOMES=0000000000001 LEN=224 CR=0\nLVACGGGSGAGDNNTPSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDYITFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPKNENNSHPGAWDVEFPDSNNPAYLRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
#         "MGYP001038205783": ">MGYP001038205783 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=1\nMKKLLALILAAALALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
#         "MGYP000111223802": ">MGYP000111223802 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=1\nMKKTLALILTAALTLSLAACGGSNGTENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001043698415": ">MGYP001043698415 PL=10 UP=0 BIOMES=0000000000001 LEN=115 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHKTLSFIINHA",
#         "MGYP001044920212": ">MGYP001044920212 PL=11 UP=0 BIOMES=0000000000001 LEN=232 CR=0\nKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYI",
#         "MGYP001045957942": ">MGYP001045957942 PL=01 UP=0 BIOMES=0000000000001 LEN=126 CR=0\nVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001051398202": ">MGYP001051398202 PL=00 UP=0 BIOMES=0000000000001 LEN=240 CR=1\nMKKLLALMLAAALALSLVACGGGGAEDNNTTPSTGNGDTTGTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCITFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPASDYTDCIMEQAYLVTDRYEYTGIPKNENNSHPGAWDVEFPDSDNPEYLRLVYFDDSIDVSQYVGKEITFSAKVMDGDYFESAFVDYIDAIIIE",
#         "MGYP001051968473": ">MGYP001051968473 PL=01 UP=0 BIOMES=0000000000001 LEN=52 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTASGGVEDSTPSMTK",
#         "MGYP001058115450": ">MGYP001058115450 PL=11 UP=0 BIOMES=0000000000001 LEN=176 CR=0\nDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSQDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNKNDSFPGAWNVEFPNAANRQYLRLVYFDDT",
#         "MGYP001061453680": ">MGYP001061453680 PL=01 UP=0 BIOMES=0000000000001 LEN=219 CR=0\nACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
#         "MGYP001064539059": ">MGYP001064539059 PL=10 UP=0 BIOMES=0000000000001 LEN=202 CR=0\nMKKLLALMLAAALALSLVACGGGNGAGDTNTSSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDD",
#         "MGYP001065068340": ">MGYP001065068340 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALVLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001067122102": ">MGYP001067122102 PL=01 UP=0 BIOMES=0000000000001 LEN=219 CR=0\nLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSIDVGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001067638354": ">MGYP001067638354 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001069194978": ">MGYP001069194978 PL=10 UP=0 BIOMES=0000000000001 LEN=224 CR=0\nALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVIVYLPSEDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNENDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
#         "MGYP001071076496": ">MGYP001071076496 PL=10 UP=0 BIOMES=0000000000001 LEN=55 CR=1\nMKKLLALILAAALALSLVACGGGGGVGDNNTPSTGNGDTTGTDTPSGGEDSTQTA",
#         "MGYP001075906476": ">MGYP001075906476 PL=00 UP=0 BIOMES=0000000000001 LEN=240 CR=1\nMKKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNNPEYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
#         "MGYP001075906285": ">MGYP001075906285 PL=11 UP=0 BIOMES=0000000000001 LEN=176 CR=1\nSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNYRRLVYFDDSIDVSQYVGKEITFSAKA",
#         "MGYP001081128765": ">MGYP001081128765 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001084165602": ">MGYP001084165602 PL=10 UP=0 BIOMES=0000000000001 LEN=62 CR=1\nMKKLLALILAAALALSLVACGGDSGAGDTNTPSGGNGDTTSTDTPSIPSGQYLNYSRIESGQ",
#         "MGYP001086225311": ">MGYP001086225311 PL=10 UP=0 BIOMES=0000000000001 LEN=74 CR=0\nMKKLLALMLAAALALSLVACGGGSGAGDNNTPSTGNGDTSSTDTPSMTKEEMLETAEDGNISELNHLIAENILS",
#         "MGYP000355976475": ">MGYP000355976475 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001108197041": ">MGYP001108197041 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=0\nMKKLLALILAAALALSLVACGGGSGAGDTNTPSGGNGDTTSTDTPSMTKEEILETAEDGDIYELNHLISENILNAKQSYCGKNLILKSTVDEIKEDCIVFNHGMGSVILYLPSQDIVNLKTNQTITIVGITNEEFVTTNESRDGGPAFDYTDCVMEQAYLVTDRYEYTGIPKNKNDSFPGAWNVEFPNAANRQYLRLVYFDDTIDVSKYEGEEIKFSAKAFRSSSSLFDYYDAIIIE",
#         "MGYP001113077337": ">MGYP001113077337 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPMAEDYTDCIMEQAYLVADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDAGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001114448909": ">MGYP001114448909 PL=00 UP=0 BIOMES=0000000000001 LEN=235 CR=0\nMKKTLALILTAALTLSLAACGGSNGTENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITTKESRYGPIAEDYTDCIMEQAYLIADRYEYSGTPKSKNSDFDGAWNVEFSNVSNPQYLRLVYFDDSVDVGQYEGKEITFSAKVIDGKYYDVVITG",
#         "MGYP001125457394": ">MGYP001125457394 PL=00 UP=0 BIOMES=0000000000001 LEN=237 CR=0\nMKKLLALMLAAALALSLVACGGGNGAGDTNTSSTGNGDTTSTDTPSMTKEEMLETAEDGNISELNHLIAENILSAKQAYCGKVLTFRGNIYGIKEDCIIFNHGKGSIIVYLPAEDIVNLKTDQWVTIIGLTNDEFITTKESRGGEPTSDYIDCIMEQAYLVTDRYEYTGIPENENNSYPGAWNVEFPDSNYRRLVYFDDSIDVSQYVGKEITFSAKAIGGDYFESAFVDYIDAIIIE",
#         "MGYP001130951151": ">MGYP001130951151 PL=01 UP=0 BIOMES=0000000000001 LEN=141 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTPSMTKEEMLEVAVNGDISELNSLIADNILSAKQAYCGKVLTMKHAVDEIKEDCIIFNHGMGSVIVYLPEEDIINLKAKQYITIVGITSEDFITT",
#         "MGYP001143804030": ">MGYP001143804030 PL=01 UP=0 BIOMES=0000000000001 LEN=48 CR=0\nMKKTLALILIAALTLSLAACGGSNGKENNNTSSTGNTTSSGGVEDSTP"
#     }
#     # Debug
#     print('There are {:d} input sequences'.format(len(fasta_sequences)))
#
#     # Make new instance of MobiDB Lite
#     mobidb = MobiDbLite(env=env)
#     # Run MobiDB Lite predictor, get disordered regions
#     disorder = mobidb.run(fasta_sequences, verbose=True)
#     # Debug
#     print('Predictions:')
#     # print(json.dumps(disorder, indent=1))
#     for acc in disorder.keys():
#         print('  Sequence accession number:', acc)
#         for i in range(len(disorder[acc])):
#             print('    [' + ', '.join([str(r) for r in disorder[acc][i]]) + ']')
#
#     # Compute threshold version of predicted disorder
#     disorder = mobidb.apply_threshold(disorder, threshold=1, inclusive=1)
#     # Debug
#     print('Thresholded predictions:')
#     for acc in disorder.keys():
#         print('  Sequence accession number:', acc)
#         for i in range(len(disorder[acc])):
#             print('    [' + ', '.join([str(r) for r in disorder[acc][i]]) + ']')
#
#     # Compute compositional bias
#     comp_bias = mobidb.compute_bias(disorder)
#     # Debug
#     print('Compositional bias: {:.02f}'.format(comp_bias))
