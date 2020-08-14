class HMMAlign(HMMER):

    # Constructor
    def __init__(self, cmd=['hmmalign'], env=os.environ.copy()):
        # Call parent constructor
        super().__init__(cmd=cmd, env=env)

    # Run HMM align
    def run(self, hmm_path, fasta_path, out_path):
        """Align an HMM model with a FASTA dataset

        Args
        hmm_path (str)          Path to HMM model file
        fasta_path (str)        Path to FASTA target sequences file
        out_path (str)          Path where to store output alignment file

        Return
        (CompletedProcess)      Result object representing a successful run

        Raise
        (CalledProcessError)    Error related to command line input
        (Error)                 Unknown error
        """
        # Define command
        cmd = self.cmd
        # Set mandatory command line options
        cmd += ['-o', out_path]  # Output file path
        cmd += [hmm_path]  # HMM model file path
        cmd += [fasta_path]   # Target dataset file path

        # Execute script
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )
