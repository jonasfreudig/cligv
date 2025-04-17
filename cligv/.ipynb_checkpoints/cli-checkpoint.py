#!/usr/bin/env python3

"""
Command-line interface for the clIGV (command line Interactive Genome Viewer).

Entry point for the command-line tool, parsing
arguments and initializing the genome browser.
"""

import argparse
import logging
import os
import sys
from typing import Optional

from rich.console import Console

from . import __version__
from .browser import GenomeBrowser


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="clIGV - command line Interactive Genome Viewer. "
        "Displays FASTA sequence, VCF variants, and BAM alignments.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "fasta", 
        help="Path to the reference genome FASTA file (indexed with .fai)"
    )
    
    # Optional arguments
    parser.add_argument(
        "-v", "--vcf", 
        help="Path to VCF file (indexed with tabix)", 
        required=False
    )
    
    parser.add_argument(
        "-b", "--bam", 
        help="Path to BAM file (indexed with .bai or .csi)", 
        required=False
    )
    
    parser.add_argument(
        "-r", "--region", 
        help="Initial region (e.g., chr1:1000-2000, chrX:50000, chrY)", 
        required=False
    )
    
    parser.add_argument(
        "-t", "--tag", 
        help="BAM tag name to use for coloring reads (e.g., 'ha' to use haplotype tags of STAR diploid)", 
        required=False
    )
    
    parser.add_argument(
        "--light-mode", 
        action='store_true',
        help="Use light color theme (for light background terminals)"
    )
    
    parser.add_argument(
        '--log-level', 
        default='INFO', 
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
        help='Set the logging level.'
    )
    
    parser.add_argument(
        '--version', 
        action='version', 
        version=f'clIGV v{__version__}'
    )
    
    return parser.parse_args()


def setup_logging(log_level: str) -> None:
    """Configure logging with the specified level."""
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - [%(funcName)s:%(lineno)d] - %(message)s',
        filename='cligv_debug.log',
        filemode='w'
    )
    logging.getLogger("pysam").setLevel(logging.WARNING)


def validate_files(fasta_path: str, vcf_path: Optional[str] = None, 
                  bam_path: Optional[str] = None) -> bool:
    """
    Validate that input files exist before attempting to open them.
    
    Returns:
        bool: True if all specified files exist, False otherwise.
    """
    console = Console()
    
    # Check FASTA file
    if not os.path.exists(fasta_path):
        console.print(f"[bold red]Error:[/bold red] FASTA file not found: {fasta_path}")
        return False
        
    # Check FASTA index
    fasta_index = fasta_path + ".fai"
    if not os.path.exists(fasta_index):
        console.print(
            f"[bold red]Error:[/bold red] FASTA index file not found: {fasta_index}.\n"
            f"Please index with 'samtools faidx {fasta_path}'"
        )
        return False
        
    # Check VCF file if specified
    if vcf_path and not os.path.exists(vcf_path):
        console.print(f"[bold red]Error:[/bold red] VCF file not found: {vcf_path}")
        return False
        
    # Check BAM file if specified
    if bam_path and not os.path.exists(bam_path):
        console.print(f"[bold red]Error:[/bold red] BAM file not found: {bam_path}")
        return False
        
    return True


def main():
    """Main entry point for the command-line tool."""
    try:
        # Parse arguments
        args = parse_args()
        
        # Check if mandatory FASTA argument is present
        if not args.fasta:
            print("Error: FASTA file argument is required.", file=sys.stderr)
            sys.exit(1)
            
        # Setup logging
        setup_logging(args.log_level)
        
        logging.info(f"Starting CLIGV v{__version__}...")
        
        # Validate input files
        if not validate_files(args.fasta, args.vcf, args.bam):
            sys.exit(1)
            
        # Initialize and run browser
        browser = GenomeBrowser(
            fasta_path=args.fasta,
            vcf_path=args.vcf,
            bam_path=args.bam,
            initial_region=args.region,
            tag_name=args.tag,
            light_mode=args.light_mode
        )
        
        browser.run()
        
    except KeyboardInterrupt:
        logging.info("Operation cancelled by user.")
        print("\nOperation cancelled by user.")
        sys.exit(0)
    except Exception as e:
        logging.exception("Unexpected error")
        print(f"\nUnexpected error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        logging.info("Browser closed.")
        print("Browser closed.")


if __name__ == "__main__":
    main()