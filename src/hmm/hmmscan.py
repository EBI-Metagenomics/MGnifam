# Dependencies
from src.hmm.hmmer import HMMER
import subprocess
import argparse
import os


class HMMScan(HMMER):

    # Constructor
    def __init__(self, cmd=['hmmsearch'], env=os.environ):
        # Call parent constructor
        super().__init__(cmd=cmd, env=env)

    # Run script
    def run(
        self, hmmdb_path='', fasta_path='',
        out_path='', tblout_path='', domtblout_path='', pfamtblout_path='',
        use_accession=False, no_alignment=False, text_width=0,
        sequence_eval=(None, None), sequence_bits=(None, None),
        domain_eval=(None, None), domain_bits=(None, None),
        cut_tc=False, cut_nc=False, cut_ga=False,
        num_sequences=None, num_domains=None, num_cpus=None
    ):
        """ Run hmmsearch

        Args
        hmmdb_path (str)        Path to HMM database path (comes from hmmpress)
        fasta_path (str)        Path to fasta sequences file
        out_path (str)          Path to output file, avoid stdout
        tblout_path (str)       Path to tblout file
        domtblout_path (str)    Path to domtblout file
        pfamtblout_path (str)   Path to pfamtblout file
        use_accession (bool)    Use accession instead of names in main output,
                                if available
        no_alignment (bool)     Avoid generating alignment in main output file,
                                greatly improves performance
        text_width (int)        Unlimit (if not set) or constraint (otherwise)
                                width of output text
        sequence_eval (tuple)   Tuple defininig reporting and inclusion
                                threshold for sequences E-value
        sequence_bits (tuple)   Tuple defining reporting and inclusion
                                threshold for sequences bit-score
        domain_eval (tuple)     Tuple defininig reporting and inclusion
                                threshold for domains E-value
        domain_bits (tuple)     Tuple defininig reporting and inclusion
                                threshold for domain bit-score
        cut_tc (bool)           Whether to use model's TC cutoff
        cut_nc (bool)           Whether to use model's NC cutoff
        cut_ga (bool)           Whether to use model's GA cutoff
        num_sequences (int)     Number of sequences compared (sequence Z-score)
        num_domains (int)       Number of domains compared (domain Z-score)
        num_cpus (int)          Number of parallel threads to use

        Return
        (CompletedProcess)      Script output

        Raise
        (CalledProcessError)    In case the script cannot be run
        (FileNotFoundError)     In case HMM database file input fasta file have
                                not been found
        """
        # Initialize command
        cmd = self.cmd
        # Define output path
        cmd += ['-o', out_path] if out_path else []
        # Define tblout path
        cmd += ['--tblout', tblout_path] if tblout_path else []
        # Define domtblout path
        cmd += ['--domtblout', domtblout_path] if domtblout_path else []
        # Define pfamtblout
        cmd += ['--pfamtblout', pfamtblout_path] if pfamtblout_path else []
        # Set use accession flag
        cmd += ['--acc'] if use_accession else []
        # Set no alignment flag
        cmd += ['--noali'] if no_alignment else []
        # Set output text width (or remove limit)
        cmd += ['--textw', text_width] if text_width else ['--notextw']
        # Set reporting and inclusion E-value thresholds for sequences
        cmd += ['-E', sequence_eval[0]] if sequence_eval[0] else []
        cmd += ['--incE', sequence_eval[1]] if sequence_eval[1] else []
        # Set reporting and inclusion E-value thresholds for domains
        cmd += ['--domE', domain_eval[0]] if domain_eval[0] else []
        cmd += ['--incdomE', domain_eval[1]] if domain_eval[1] else []
        # Set reporting and inclusion bit-score threshold for sequences
        cmd += ['-T', sequence_bits[0]] if sequence_bits[0] else []
        cmd += ['--incT', sequence_bits[1]] if sequence_bits[1] else []
        # Set reporting and inclusion bit-score threshold for domains
        cmd += ['--domT', domain_bits[0]] if domain_bits[0] else []
        cmd += ['--incdomT', domain_bits[1]] if domain_bits[1] else []
        # Set gathering thresholds
        cmd += ['--cut_ga'] if cut_ga else []
        cmd += ['--cut_nc'] if cut_nc else []
        cmd += ['--cut_tc'] if cut_tc else []
        # Set Z-scores for sequences and domains
        cmd += ['-Z', num_sequences] if num_sequences else []
        cmd += ['--domZ', num_domains] if num_domains else []
        # Set number of CPUs (threads)
        cmd += ['--cpu', num_cpus] if num_cpus else []
        # Add mandatory parameters: HMM database and FASTA file
        cmd += [hmmdb_path, fasta_path]
        # Ensure that each parameter is string
        cmd = [str(cmd[i]) for i in range(len(cmd))]

        # Run hmmscan, return CompletedProcess instance
        return subprocess.run(
            capture_output=True,  # Capture console output
            check=True,  # Check that script ran successfully
            encoding='utf-8',  # Set output encoding
            env=self.env,  # Set script environment
            args=cmd  # Set command line arguments
        )


# Run hmmscan from command line
if __name__ == '__main__':

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Run hmmscan')
    # Initialize hmmscan command
    parser.add_argument(
        '--cmd', nargs='+', type=str, default=['hmmscan'],
        help='Path to hmmscan executable'
    )
    # Define path to input HMM database (directory)
    parser.add_argument(
        '--hmmdb_path', type=str, required=True,
        help='Path to HMM database directory (requires hmmpress)'
    )
    # Define path to input FASTA file
    parser.add_argument(
        '--fasta_path', type=str, required=True,
        help='Path to input FASTA (consensus) file'
    )
    # Define output path
    parser.add_argument(
        '-o', '--out_path', type=str, default='',
        help='Output file path'
    )
    # Define tblout path
    parser.add_argument(
        '--tblout_path', type=str, default='',
        help='Path to tblout output file, if any'
    )
    # Define domtblout path
    parser.add_argument(
        '--domtblout_path', type=str, default='',
        help='Path to domtblout output file, if any'
    )
    # Define pfamtblout path
    parser.add_argument(
        '--pfamtblout_path', type=str, default='',
        help='Path to pfamtblout output file, if any'
    )
    # Define accession flag
    parser.add_argument(
        '--use_accession', type=int, default=0,
        help='Use accession numbers instead of names'
    )
    # Define noalign flag
    parser.add_argument(
        '--no_alignment', type=int, default=1,
        help='Make sequences alignments (at a very high computational cost)'
    )
    # Define output text width
    parser.add_argument(
        '--text_width', type=int, default=0,
        help='Maximum output text width, unlimited by default'
    )
    # Define e-value reporting threshold for sequences
    parser.add_argument(
        '--reporting_seqe', type=float, required=False,
        help='Sequence e-value reporting threshold'
    )
    # Define e-value inclusion threshold for sequences
    parser.add_argument(
        '--inclusion_seqe', type=float, required=False,
        help='Sequence e-value inclusion threshold'
    )
    # Define e-value reporting threshold for domains
    parser.add_argument(
        '--reporting_dome', type=float, required=False,
        help='Domain e-value reporting threshold'
    )
    # Define e-value inclusion threshold for domains
    parser.add_argument(
        '--inclusion_dome', type=float, required=False,
        help='Domain e-value inclusion threshold'
    )
    # Define bit-score reporting threshold for sequences
    parser.add_argument(
        '--reporting_seqb', type=float, required=False,
        help='Sequence bit-score reporting threshold'
    )
    # Define bit-score inclusion threshold for sequences
    parser.add_argument(
        '--inclusion_seqb', type=float, required=False,
        help='Sequence bit-score inclusion threshold'
    )
    # Define bit-score reporting threshold for domains
    parser.add_argument(
        '--reporting_domb', type=float, required=False,
        help='Domain bit-score reporting threshold'
    )
    # Define bit-score inclusion threshold for domains
    parser.add_argument(
        '--inclusion_domb', type=float, required=False,
        help='Domain bit-score inclusion threshold'
    )
    # Define gathering threshold (GA) flag
    parser.add_argument(
        '--cut_ga', type=bool, required=False,
        help='Whether to use model GA cutoff'
    )
    # Define TC flag
    parser.add_argument(
        '--cut_tc', type=bool, required=False,
        help='Whether to use model TC cutoff'
    )
    # Define NC flag
    parser.add_argument(
        '--cut_nc', type=bool, required=False,
        help='Whether to use model NC cutoff'
    )
    # Define number of sequences
    parser.add_argument(
        '--num_sequences', type=int, default=0,
        help='Number of sequences compared (sequences Z-score)'
    )
    # Define number of domains
    parser.add_argument(
        '--num_domains', type=int, default=0,
        help='Number of domains compared (domains Z-score)'
    )
    # Define number of CPUs
    parser.add_argument(
        '--num_cpus', type=int, default=0,
        help='Number of CPUs used in parallel for multithreading'
    )

    # Retrieve arguments
    args = parser.parse_args()

    # Try running an HMMScan instance
    try:
        # Instantiate HMMScan object
        hmmscan = HMMScan(cmd=args.cmd, env=os.environ)
        # Run HMMScan instance with given arguments
        ran = hmmscan(
            hmmdb_path=args.hmmdb_path,
            fasta_path=args.fasta_path,
            out_path=args.out_path,
            tblout_path=args.tblout_path,
            domtblout_path=args.domtblout_path,
            pfamtblout_path=args.pfamtblout_path,
            use_accession=bool(args.use_accession),
            no_alignment=bool(args.no_alignment),
            text_width=args.text_width,
            sequence_eval=(args.reporting_seqe, args.inclusion_seqe),
            sequence_bits=(args.reporting_seqb, args.inclusion_seqb),
            domain_eval=(args.reporting_dome, args.inclusion_dome),
            domain_bits=(args.reporting_domb, args.inclusion_domb),
            cut_ga=args.cut_ga, cut_tc=args.cut_tc, cut_nc=args.cut_nc,
            num_sequences=args.num_sequences,
            num_domains=args.num_domains,
            num_cpus=args.num_cpus
        )

        # Case no output path has been set, return STDOUT
        if not args.out_path:
            # Verbose
            print('Successfully ran hmmscan:')
            print(ran.stdout.strip())

    # Intercept exception
    except subprocess.CalledProcessError as err:
        # Show error STDERR
        raise Exception('\n'.join([
            'Error: hmmscan returned code {:}:'.format(err.returncode),
            err.stderr.strip()
        ]))
