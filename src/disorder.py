# Dependencies
import numpy as np
import subprocess
import tempfile
import os

class MobiDBLite(object):

    # Constructor
    def __init__(self, cmd=['mobidb_lite.py'], env=os.environ.copy()):
        # Save attributes
        self.cmd = cmd  # Path to MobiDB Lite executable
        self.env = env  # Environmental variables

    # Run method
    def run(self, sequences, acc_regex='^>(\S+)', ret_comp_bias=False, verbose=True):
        """Run MobiDB Lite disorder predictor
        Takes as input a list of sequences and returns a lis of boolean (1/0)
        values, where an item has value 1 if the residue in the same position is
        predicted to be disordered, 0 otherwise.

        Args
        sequences (list)    List of strings containing fasta entries
        verbose (bool)      Wether to show verbose log or not

        Return
        (list(list))        List of boolean (1/0) lists
        """
        # Define temporary input file
        fasta_file = tempfile.NamedTemporaryFile(delete=False)
        # Define temporary output file
        out_file = tempfile.NamedTemporaryFile(delete=False)

        # Make input fasta file
        with open(fasta_file.name, 'w') as ff:
            # Loop through each input sequence
            for seq in sequences:
                # Write out fasta entry
                ff.write(seq + '\n')

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

        # Debug
        print(cmd_out)

        with open(out_file.name, 'r') as of:
            # Debug
            print(of.read())

        # Delete temporary files
        os.remove(fasta_file.name)
        os.remove(out_file.name)

    @staticmethod
    def compute_comp_bias(sequences, threshold=1):
        raise NotImplementedError

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

    # Define list of test fasta entries
    sequences = [
        ">MGYP001026457211 PL=10 UP=0 BIOMES=0000000011000 LEN=188 CR=0\nNTTCENADRLGTIPADSYAKVNSPSFLGIPLTTTPKPDAPSTQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWESVFSEVALSINALNAAVDQPTSPNDYSGQLKFTGKRKITALNLTNIKTTTSEYATVIGMRADNKELAYEFICIDNYIYMRTGKGDTWNNAISIISKFTSF",
        ">MGYP000848664103 PL=00 UP=0 BIOMES=0000000011000 LEN=297 CR=1\nMAEYSSELDKITYAELALSLQNTIKNNLAHTKDQVIHVTQEDKNKWNQISDIPEATETKKGALTPQEKIKLKNIEERANNYTHPTSGVTAGQYIQVEVNAEGHVVAGHNPTKINTTCENADRLGTIPADSYAKVNSPSFLGIPLTTTPKPDAPSTQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWQSVFSEVALSINALNAAVDQPTSPNDYSGQLKFTGKRKITALKLTNIKATTSEYATVIGMRADNRELAYEFICIDNYIYMRTGKGDTWNNAISIIKD",
        ">MGYP001028524128 PL=01 UP=0 BIOMES=0000000011000 LEN=26 CR=0\nMAEYSSELDKITYAELALSLQNTIKN",
        ">MGYP000854299633 PL=00 UP=0 BIOMES=0000000011000 LEN=297 CR=0\nMAEYSSELDKITYAELALSLQNTIKSNLAHTKDQVIHVTQEDKNKWNQISDIPEATETKKGALTPQEKIKLKNIEERANNYTHPTSGVTAGQYIQVEVNAEGHVVAGHNPTKINTTCENADRLGTIPADSYAKVNSPSFLGIPLTPTPKLDAPASQIVNIEYLNSQPTYIRQKTAPEKALSGKLWIGNNNCLNAYNNDGWQSVFSEVALSINALNTAVDQPTSPNDYSGQFKFTGKRKITALKLTNIKATTSEYATVIGMRADNRELAYEFICIDNYIYMRTGKGDTWNNAISIIKD"
    ]
    # Make new instance of MobiDB Lite
    mdb_lite = MobiDBLite(env=env)
    # Run MobiDB Lite predictor
    mdb_lite.run(sequences)
