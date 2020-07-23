# Dependencies
import argparse
import sys
import os

# Set path for custom dependencies
sys.path.append(os.path.dirname(__file__) + '/..')

# Custom dependencies
from src.pipeline.release import ReleasePipeline


# Main
if __name__ == '__main__':

    # Define new argument parser
    parser = argparse.ArgumentParser(description='Make new MGnifam release')
    # Feed command line arguments to it
    parser.add_argument(
        '-c', '--cores', type=int, default=1,
        help='Maximum number of jobs per core'
    )
    parser.add_argument(
        '-m', '--memory', type=str, default='8 GB',
        help='Maximum amount of memory per job'
    )
    parser.add_argument(
        '--cluster', type=str, default='LSFCluster',
        help='Type of cluster to use (only LSFCluster implemented yet)'
    )
    # Parse the arguments defined above
    args = parser.parse_args()

    
