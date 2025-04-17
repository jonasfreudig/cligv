"""
Main browser module for the command line Interactive Genome Viewer (clIGV).

This module contains the GenomeBrowser class that renders and manages the
interactive genome viewer interface.
"""

import logging
import math
import os
import sys
import time
from typing import Dict, List, Optional, Set, Tuple, Any

# Core libraries
try:
    import pysam
except ImportError:
    print("Error: pysam library not found. Please install it: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    from rich.console import Console
    from rich.layout import Layout
    from rich.panel import Panel
    from rich.text import Text
    from rich.live import Live
    from rich import box
except ImportError:
    print("Error: rich library not found. Please install it: pip install rich", file=sys.stderr)
    sys.exit(1)

from .alignment import AlignmentRead
from .themes import get_theme, get_color
from .utils import (
    DEFAULT_REGION_WIDTH,
    MAX_REGION_WIDTH,
    MAX_FETCH_SIZE,
    format_position, 
    genomic_to_screen_coord,
    parse_region
)

# --- Constants ---
MAX_ALIGNMENT_ROWS = 25  # Max rows to display
ZOOM_FACTOR = 1.5  # Zoom factor
SEQUENCE_DISPLAY_THRESHOLD = 200  # Switch to condensed view above this width
DETAILED_VIEW_THRESHOLD = 1000  # Show detailed variant info below this width
VARIANT_DISPLAY_THRESHOLD = 1000  # Show variant alleles below this width
NAVIGATION_OVERLAP = 0.5  # Navigation overlap (50%)
COVERAGE_HEIGHT = 10  # Height of coverage track


class GenomeBrowser:
    """
    Manages the state and rendering of the command-line genome browser with proper scaling.
    """
    
    def __init__(self, fasta_path: str, vcf_path: Optional[str] = None, 
                 bam_path: Optional[str] = None, initial_region: Optional[str] = None, 
                 tag_name: Optional[str] = None, light_mode: bool = False):
        """
        Initialize the genome browser with files and initial region.
        
        Args:
            fasta_path: Path to the reference FASTA file
            vcf_path: Optional path to a VCF file
            bam_path: Optional path to a BAM file
            initial_region: Optional initial region to display
            tag_name: Optional BAM tag name to use for coloring reads
            light_mode: Use light color theme if True
        """
        self.console = Console()
        logging.info("Initializing GenomeBrowser...")
        self.fasta_file = None
        self.vcf_file = None
        self.bam_file = None
        self.tag_name = tag_name
        self.light_mode = light_mode
        
        # Set up theme
        self.theme = get_theme(light_mode=self.light_mode)
        logging.info(f"Using {'light' if light_mode else 'dark'} theme")
        
        # Open files
        self._open_files(fasta_path, vcf_path, bam_path)
        
        # Initialize state
        self.chrom = None
        self.start = None
        self.end = None
        self.current_width = self.console.width
        self.reference_sequence = ""
        self.current_variants = {}  # Position -> variant info
        self.processed_alignments = []  # Processed alignment objects
        
        # Create layout
        self.layout = self._create_layout()
        logging.debug("Layout created.")
        
        # Set initial region
        self._set_initial_region(initial_region)
        
    def _open_files(self, fasta_path, vcf_path, bam_path):
        """
        Open genome files with proper error handling.
        
        Args:
            fasta_path: Path to reference FASTA file
            vcf_path: Optional path to VCF file
            bam_path: Optional path to BAM file
        """
        # Open FASTA (required)
        try:
            fasta_index = fasta_path + ".fai"
            if not os.path.exists(fasta_index):
                raise FileNotFoundError(f"FASTA index file not found: {fasta_index}. "
                                       f"Please index with 'samtools faidx'.")
            self.fasta_file = pysam.FastaFile(fasta_path)
            logging.info(f"Opened FASTA: {fasta_path}")
            self.console.print(f"[{get_color(self.theme, 'SUCCESS_STYLE')}]Opened FASTA:[/] {fasta_path}")
        except Exception as e:
            logging.exception(f"Error opening FASTA file '{fasta_path}'")
            self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error opening FASTA file '{fasta_path}':[/] {e}")
            sys.exit(1)

        # Open VCF (optional)
        if vcf_path:
            try:
                # Find VCF index
                vcf_index = vcf_path + ".tbi"
                if not os.path.exists(vcf_index):
                    if vcf_path.endswith(".vcf.gz"):
                        alt_vcf_index = vcf_path + ".tbi"
                    elif vcf_path.endswith(".vcf"):
                        alt_vcf_index = vcf_path + ".gz.tbi"
                    else:
                        alt_vcf_index = vcf_path + ".vcf.gz.tbi"
                    
                    if os.path.exists(alt_vcf_index):
                        vcf_index = alt_vcf_index
                    else:
                        raise FileNotFoundError(f"VCF index file not found. "
                                               f"Please index with 'tabix -p vcf {vcf_path}'.")
                
                self.vcf_file = pysam.VariantFile(vcf_path)
                logging.info(f"Opened VCF: {vcf_path}")
                self.console.print(f"[{get_color(self.theme, 'SUCCESS_STYLE')}]Opened VCF:[/] {vcf_path}")
            except Exception as e:
                logging.warning(f"Warning opening VCF file '{vcf_path}': {e}", exc_info=False)
                self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Warning opening VCF file '{vcf_path}':[/] {e}")
                self.vcf_file = None

        # Open BAM (optional)
        if bam_path:
            try:
                # Find BAM index
                bai_path = bam_path + ".bai"
                csi_path = bam_path + ".csi"
                
                if os.path.exists(bai_path):
                    index_path = bai_path
                elif os.path.exists(csi_path):
                    index_path = csi_path
                else:
                    raise FileNotFoundError(f"BAM index file not found. "
                                          f"Please index with 'samtools index {bam_path}'.")
                
                self.bam_file = pysam.AlignmentFile(bam_path, "rb")
                logging.info(f"Opened BAM: {bam_path}")
                self.console.print(f"[{get_color(self.theme, 'SUCCESS_STYLE')}]Opened BAM:[/] {bam_path}")
                
                if self.tag_name:
                    self.console.print(f"[{get_color(self.theme, 'SUCCESS_STYLE')}]Will color reads by tag:[/] {self.tag_name}")
            except Exception as e:
                logging.warning(f"Warning opening BAM file '{bam_path}': {e}", exc_info=False)
                self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Warning opening BAM file '{bam_path}':[/] {e}")
                self.bam_file = None
    
    def _create_layout(self) -> Layout:
        """
        Create the Rich layout structure with all necessary tracks.
        
        Returns:
            Layout: Configured Rich layout
        """
        layout = Layout()
        
        # Define the main top-level splits
        layout.split(
            Layout(name="header", size=3),
            Layout(name="main", ratio=1),
            Layout(name="footer", size=3)
        )
        
        # Split the main area for different tracks
        layout["main"].split_column(
            Layout(name="overview", size=3),
            Layout(name="ruler", size=1),
            Layout(name="sequence", size=3),
            Layout(name="variants", size=5, visible=bool(self.vcf_file)),
            Layout(name="coverage", size=COVERAGE_HEIGHT, visible=bool(self.bam_file)),
            Layout(name="alignments", ratio=1, visible=bool(self.bam_file)),
        )
        
        return layout
    
    def _set_initial_region(self, initial_region):
        """
        Set the initial viewing region.
        
        Args:
            initial_region: Region string to set
        """
        logging.info(f"Setting initial region from input: {initial_region}")
        parsed_region = parse_region(initial_region)
        
        if parsed_region:
            chrom, start, end = parsed_region
            if not self.set_region(chrom, start, end):
                logging.error("Initial region is invalid or chromosome not found.")
                self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error:[/] "
                                 f"Initial region is invalid or chromosome not found.")
                sys.exit(1)
        elif self.fasta_file and self.fasta_file.references:
            first_chrom = self.fasta_file.references[0]
            logging.info(f"No initial region provided, using first reference: {first_chrom}")
            if not self.set_region(first_chrom, 1, DEFAULT_REGION_WIDTH):
                logging.error("Could not set default initial region.")
                self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error:[/] "
                                 f"Could not set default initial region.")
                sys.exit(1)
        else:
            logging.error("No initial region provided and FASTA file seems empty or invalid.")
            self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error:[/] "
                             f"No initial region provided and FASTA file seems empty or invalid.")
            sys.exit(1)
            
        logging.info(f"GenomeBrowser initialized. Region: {self.chrom}:{self.start}-{self.end}")
    
    def set_region(self, chrom: str, start: int, end: int) -> bool:
        """
        Set the current genomic region being viewed.
        Validates coordinates and fetches reference sequence.
        
        Args:
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based)
            
        Returns:
            bool: True on success, False on failure
        """
        logging.info(f"Attempting to set region: {chrom}:{start}-{end}")
        
        if not self.fasta_file:
            logging.error("FASTA file not loaded in set_region.")
            self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error:[/] FASTA file not loaded.")
            return False

        # Check if chromosome exists
        if chrom not in self.fasta_file.references:
            logging.warning(f"Chromosome '{chrom}' not found in FASTA file.")
            self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Warning:[/] "
                             f"Chromosome '{chrom}' not found in FASTA file. "
                             f"Available: {', '.join(self.fasta_file.references[:5])}...")
            return False

        # Validate coordinates
        chrom_len = self.fasta_file.get_reference_length(chrom)
        new_start = max(1, start)  # 1-based start
        new_end = min(chrom_len, end)

        if new_start > new_end:
            logging.warning(f"Invalid region coordinates: {new_start}-{new_end}")
            self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Warning:[/] "
                             f"Invalid region coordinates after clamping. Keeping previous region.")
            return False

        # Update state
        self.chrom = chrom
        self.start = new_start
        self.end = new_end
        logging.info(f"Region set to: {self.chrom}:{self.start}-{self.end}")

        # Fetch reference sequence
        try:
            fetch_start = self.start - 1  # 0-based for pysam
            fetch_end = self.end  # exclusive end

            if fetch_start < 0:
                fetch_start = 0

            if fetch_start >= fetch_end:
                if fetch_start < chrom_len:
                    fetch_end = fetch_start + 1
                else:
                    raise ValueError(f"Region start {fetch_start} beyond chromosome length {chrom_len}")

            logging.debug(f"Fetching sequence: {self.chrom}:{fetch_start}-{fetch_end}")
            self.reference_sequence = self.fasta_file.fetch(self.chrom, fetch_start, fetch_end)

            # Validate sequence length
            expected_len = self.end - self.start + 1
            if len(self.reference_sequence) < expected_len:
                self.reference_sequence = self.reference_sequence.ljust(expected_len, 'N')
            elif len(self.reference_sequence) > expected_len:
                self.reference_sequence = self.reference_sequence[:expected_len]

            # Clear caches when region changes
            self.current_variants = {}
            self.processed_alignments = []

        except ValueError as e:
            logging.error(f"ValueError fetching sequence: {e}")
            self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error fetching sequence: {e}[/]")
            self.reference_sequence = "?" * (self.end - self.start + 1)
            return False
        except Exception as e:
            logging.exception(f"Unexpected error fetching sequence")
            self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Unexpected error fetching sequence: {e}[/]")
            self.reference_sequence = "?" * (self.end - self.start + 1)
            return False

        return True
    
    def zoom(self, zoom_in=True) -> bool:
        """
        Zoom in or out on the current region.
        
        Args:
            zoom_in: If True, zoom in, otherwise zoom out
            
        Returns:
            bool: True on success, False on failure
        """
        if not self.chrom or not self.start or not self.end:
            return False

        region_width = self.end - self.start + 1
        center = self.start + region_width // 2

        # Calculate new width
        if zoom_in:
            new_width = max(10, int(region_width / ZOOM_FACTOR))
        else:
            new_width = int(region_width * ZOOM_FACTOR)
            
            # Check if width is getting too large
            if new_width > MAX_REGION_WIDTH:
                self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Warning:[/] "
                                 f"Approaching large region size ({format_position(new_width)} bp).")
                if new_width > MAX_REGION_WIDTH * 2:
                    new_width = MAX_REGION_WIDTH
                    self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Limited zoom to "
                                     f"{format_position(new_width)} bp for performance.[/]")

        # Calculate new start/end
        new_start = max(1, center - new_width // 2)
        new_end = new_start + new_width - 1

        logging.info(f"Zooming {'in' if zoom_in else 'out'}: {self.chrom}:{self.start}-{self.end} -> "
                    f"{self.chrom}:{new_start}-{new_end}")
        return self.set_region(self.chrom, new_start, new_end)
        
    def navigate(self, forward=True) -> bool:
        """
        Navigate forward or backward with overlap.
        
        Args:
            forward: If True, move forward, otherwise move backward
            
        Returns:
            bool: True on success, False on failure
        """
        if not self.chrom or not self.start or not self.end:
            return False
            
        # Calculate move amount with overlap
        region_width = self.end - self.start + 1
        move_amount = int(region_width * (1 - NAVIGATION_OVERLAP))
        
        if not forward:
            move_amount = -move_amount
            
        # Calculate new positions
        new_start = max(1, self.start + move_amount)
        new_end = new_start + region_width - 1
        
        logging.info(f"Navigating {'forward' if forward else 'backward'}: "
                    f"{self.chrom}:{self.start}-{self.end} -> {self.chrom}:{new_start}-{new_end}")
        return self.set_region(self.chrom, new_start, new_end)
    
    def toggle_theme(self) -> None:
        """Toggle between light and dark themes."""
        self.light_mode = not self.light_mode
        self.theme = get_theme(light_mode=self.light_mode)
        logging.info(f"Switched to {'light' if self.light_mode else 'dark'} theme")
        
        # Clear processed data that depends on theme colors
        self.processed_alignments = []
        
    def _process_alignments(self) -> List[AlignmentRead]:
        """
        Process BAM alignments into AlignmentRead objects with proper variant highlighting.
        
        Returns:
            List[AlignmentRead]: Processed alignment objects
        """
        if not self.bam_file:
            return []
        
        if self.processed_alignments:
            return self.processed_alignments  # Use cached results if available
            
        logging.debug(f"Processing alignments for {self.chrom}:{self.start}-{self.end}")
        
        view_width = self.end - self.start + 1
        view_start_0based = self.start - 1
        ref_seq_upper = self.reference_sequence.upper()
        
        # Define tag colors with theme-specific values
        tag_colors = {
            "0": "bold black on green",
            "1": "bold black on sky_blue1",
            "2": "bold black on dark_orange"
        }
        default_tag_color = "bold white on grey30"
        
        processed_reads = []
        
        try:
            # Define fetch coordinates
            fetch_start = self.start - 1
            fetch_end = self.end
            
            if fetch_start < 0:
                fetch_start = 0
                
            # Skip fetch if region is invalid or too large
            if fetch_start >= fetch_end:
                logging.debug("Fetch region has zero width.")
                return []
                
            if self.end - self.start > MAX_FETCH_SIZE:
                logging.warning(f"Region exceeds max fetch size ({MAX_FETCH_SIZE}bp)")
                return []
                
            # Fetch reads
            logging.debug(f"Fetching BAM records")
            try:
                fetched_reads = self.bam_file.fetch(self.chrom, fetch_start, fetch_end)
            except Exception as e:
                logging.error(f"Error fetching BAM records: {e}")
                return []
                
            # Process reads
            for read in fetched_reads:
                # Skip secondary alignments
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                    
                # Get tag value if specified
                tag_value = None
                tag_style = None
                if self.tag_name:
                    try:
                        if read.has_tag(self.tag_name):
                            tag_value = read.get_tag(self.tag_name)
                            tag_str = str(tag_value)
                            tag_style = tag_colors.get(tag_str, default_tag_color)
                    except Exception as e:
                        logging.warning(f"Error getting tag: {e}")
                
                # Create alignment object
                alignment = AlignmentRead(
                    read=read,
                    ref_start=read.reference_start,
                    ref_end=read.reference_end or read.reference_start + 1,
                    has_tag=tag_value is not None,
                    tag_value=tag_value,
                    tag_style=tag_style
                )
                
                try:
                    q_seq = read.query_sequence
                    ref_idx = read.reference_start
                    q_idx = 0
                    
                    # Get CIGAR operations
                    cigar_tuples = read.cigartuples
                    if not cigar_tuples:
                        continue
                        
                    # Process CIGAR operations
                    for cigar_op, cigar_len in cigar_tuples:
                        if cigar_op == 0:  # Match/Mismatch
                            for i in range(cigar_len):
                                current_ref_pos = ref_idx + i
                                
                                # Check if position is in view
                                if view_start_0based <= current_ref_pos < view_start_0based + view_width:
                                    view_idx = current_ref_pos - view_start_0based
                                    
                                    # Process match/mismatch
                                    if (q_idx + i) < len(q_seq or "") and view_idx < len(ref_seq_upper):
                                        q_base = q_seq[q_idx + i].upper() if q_seq else '?'
                                        ref_base = ref_seq_upper[view_idx]
                                        
                                        # Check variant
                                        ref_pos_1based = current_ref_pos + 1
                                        at_variant = ref_pos_1based in self.current_variants
                                        
                                        if q_base == ref_base:
                                            base_char = get_color(self.theme, "MATCH_CHAR")
                                            if tag_style and not at_variant:
                                                base_style = tag_style
                                            else:
                                                match_style = get_color(self.theme, "MATCH_STYLE_FWD" if not read.is_reverse else "MATCH_STYLE_REV")
                                                base_style = match_style
                                        else:
                                            base_char = q_base
                                            
                                            # Check if mismatch matches a variant
                                            if at_variant:
                                                var_info = self.current_variants[ref_pos_1based]
                                                if q_base in var_info['alt']:
                                                    base_style = get_color(self.theme, "VARIANT_ALLELE_STYLE")  # Variant match
                                                else:
                                                    base_style = get_color(self.theme, "MISMATCH_STYLE")
                                            else:
                                                base_style = get_color(self.theme, "MISMATCH_STYLE")
                                        
                                        # Add to alignment representation
                                        alignment.add_position(current_ref_pos, base_char, base_style)
                                    else:
                                        alignment.add_position(current_ref_pos, "?", "dim")
                            
                            # Update indices
                            ref_idx += cigar_len
                            q_idx += cigar_len
                            
                        elif cigar_op == 1:  # Insertion
                            # Mark insertion at previous position
                            prev_ref_pos = ref_idx - 1
                            if view_start_0based <= prev_ref_pos < view_start_0based + view_width:
                                insertion_char = get_color(self.theme, "INSERTION_CHAR")
                                insertion_style = get_color(self.theme, "INSERTION_STYLE")
                                alignment.add_position(prev_ref_pos, insertion_char, insertion_style)
                            
                            # Advance query index
                            q_idx += cigar_len
                            
                        elif cigar_op == 2:  # Deletion
                            # Mark deletion along reference
                            for i in range(cigar_len):
                                current_ref_pos = ref_idx + i
                                if view_start_0based <= current_ref_pos < view_start_0based + view_width:
                                    deletion_char = get_color(self.theme, "DELETION_CHAR")
                                    deletion_style = get_color(self.theme, "DELETION_STYLE")
                                    alignment.add_position(current_ref_pos, deletion_char, deletion_style)
                            
                            # Advance reference index
                            ref_idx += cigar_len
                            
                        elif cigar_op == 3:  # Skip (e.g., intron)
                            # Show skip along reference
                            for i in range(cigar_len):
                                current_ref_pos = ref_idx + i
                                if view_start_0based <= current_ref_pos < view_start_0based + view_width:
                                    alignment.add_position(current_ref_pos, " ", "dim")
                            
                            # Advance reference index
                            ref_idx += cigar_len
                            
                        elif cigar_op == 4:  # Soft clip
                            q_idx += cigar_len
                            
                        elif cigar_op == 5:  # Hard clip
                            pass
                            
                        elif cigar_op == 6:  # Padding
                            pass
                    
                    # Add read if it has any positions
                    if alignment.representation:
                        processed_reads.append(alignment)
                        
                except Exception as e:
                    logging.warning(f"Error processing read {read.query_name}: {e}")
            
            logging.debug(f"Processed {len(processed_reads)} reads")
            
        except Exception as e:
            logging.exception(f"Error in alignment processing: {e}")
        
        # Sort reads by start position for consistent display
        processed_reads.sort(key=lambda r: r.ref_start)
        self.processed_alignments = processed_reads
        return processed_reads
    
    def _calculate_coverage(self, alignments, view_start, view_width):
        """
        Calculate coverage at each genomic position.
        
        Args:
            alignments: List of processed alignments
            view_start: 0-based start position of the view
            view_width: Width of the view in bases
            
        Returns:
            list: Coverage values for each position
        """
        coverage = [0] * view_width
        
        # Count alignment coverage at each position
        for aln in alignments:
            for pos, char, _ in aln.representation:
                rel_pos = pos - view_start
                if 0 <= rel_pos < view_width:
                    coverage[rel_pos] += 1
                    
        return coverage
    
    def _scale_coverage(self, coverage, view_width, screen_width):
        """
        Scale coverage to fit screen width.
        
        Args:
            coverage: List of coverage values
            view_width: Width of the view in bases
            screen_width: Available screen width
            
        Returns:
            list: Scaled coverage values
        """
        if view_width <= screen_width:
            # No scaling needed
            return coverage[:screen_width]
            
        # Apply scaling
        scale_factor = view_width / screen_width
        scaled_coverage = []
        
        for i in range(screen_width):
            start_idx = int(i * scale_factor)
            end_idx = int((i + 1) * scale_factor)
            # Average coverage in this section
            section = coverage[start_idx:end_idx]
            if section:
                scaled_coverage.append(sum(section) / len(section))
            else:
                scaled_coverage.append(0)
                
        return scaled_coverage
    
    def _get_coverage_display(self, coverage, max_height=COVERAGE_HEIGHT):
        """
        Convert coverage data to display rows.
        
        Args:
            coverage: List of coverage values
            max_height: Maximum height of the coverage display
            
        Returns:
            tuple: (display_rows, max_coverage)
        """
        if not coverage:
            return [], 0
            
        # Find maximum coverage for scaling
        max_coverage = max(coverage) if max(coverage) > 0 else 1
        
        # Generate display lines
        display_rows = []
        for height in range(max_height, 0, -1):
            row = []
            threshold = max_coverage * height / max_height
            
            for cov in coverage:
                if cov >= threshold:
                    row.append(('█', get_color(self.theme, "COVERAGE_BAR")))
                else:
                    row.append((' ', 'default'))
                    
            display_rows.append(row)
            
        return display_rows, max_coverage
    
    def _layout_alignments(self, alignments, view_start, view_width, screen_width, max_rows):
        """
        Layout alignments into rows without overlaps, with proper scaling.
        
        Args:
            alignments: List of processed alignments
            view_start: 0-based start position of the view
            view_width: Width of the view in bases
            screen_width: Available screen width
            max_rows: Maximum number of rows to display
            
        Returns:
            list: Display rows for alignments
        """
        if not alignments:
            return []
            
        # Track the end position of each row
        row_ends = []
        read_assignments = []  # (row_idx, read) pairs
        
        # Assign reads to rows
        for aln in alignments:
            # Find a row where this read fits
            assigned = False
            for row_idx, end_pos in enumerate(row_ends):
                if aln.ref_start > end_pos:
                    # This read fits in this row
                    row_ends[row_idx] = max(end_pos, aln.ref_end)
                    read_assignments.append((row_idx, aln))
                    assigned = True
                    break
                    
            # If no existing row works, add a new row
            if not assigned and len(row_ends) < max_rows:
                row_ends.append(aln.ref_end)
                read_assignments.append((len(row_ends) - 1, aln))
                
        # Sort by row for stable display
        read_assignments.sort(key=lambda x: x[0])
        
        # Generate display rows
        display_rows = []
        for row_idx in range(min(len(row_ends), max_rows)):
            # Create empty row
            row_display = [(' ', 'default')] * screen_width
            
            # Add reads assigned to this row
            for assigned_row, aln in read_assignments:
                if assigned_row == row_idx:
                    # Generate display for this read
                    read_display = aln.generate_display(view_start, view_width, screen_width)
                    
                    # Merge into row (read_display overwrites row_display)
                    for i, (char, style) in enumerate(read_display):
                        if char != ' ':  # Only overwrite with non-space chars
                            row_display[i] = (char, style)
                            
            display_rows.append(row_display)
            
        # Add ellipsis if we have more rows than we can show
        if len(row_ends) > max_rows:
            ellipsis_row = [(' ', 'default')] * screen_width
            ellipsis_pos = min(screen_width - 3, screen_width // 2)
            ellipsis_row[ellipsis_pos] = ('.', get_color(self.theme, "WARNING_STYLE"))
            ellipsis_row[ellipsis_pos + 1] = ('.', get_color(self.theme, "WARNING_STYLE"))
            ellipsis_row[ellipsis_pos + 2] = ('.', get_color(self.theme, "WARNING_STYLE"))
            display_rows.append(ellipsis_row)
            
        return display_rows
    
    def _fetch_variants(self):
        """
        Fetch variants in the current region from VCF.
        
        Returns:
            dict: Dictionary of variant information by position
        """
        if not self.vcf_file:
            return {}
            
        variants = {}
        
        try:
            # Fetch variants in region
            for record in self.vcf_file.fetch(self.chrom, self.start - 1, self.end):
                pos = record.pos  # 1-based position
                
                # Only process if within view
                if self.start <= pos <= self.end:
                    # Format variant info
                    alts_str = ','.join(str(a) for a in record.alts) if record.alts else '.'
                    
                    variants[pos] = {
                        'ref': record.ref,
                        'alt': alts_str,
                        'position': pos
                    }
                    
                    # Add additional details if available
                    if hasattr(record, 'qual') and record.qual is not None:
                        variants[pos]['qual'] = record.qual
                        
                    # Add important INFO fields
                    for key, value in record.info.items():
                        if key in ['DP', 'AF', 'AC']:
                            variants[pos][key] = value
                            
        except Exception as e:
            logging.warning(f"Error fetching variants: {e}")
            
        self.current_variants = variants
        return variants
    
    def _get_sequence_track(self) -> Panel:
        """
        Format the sequence track with proper scaling.
        
        Returns:
            Panel: Formatted sequence panel
        """
        logging.debug(f"Getting sequence track. Ref seq len: {len(self.reference_sequence)}")
        
        display_width = self.end - self.start + 1
        
        # Calculate exact available width for text within the panel
        panel_borders_and_padding = 4  # 2 for borders, 2 for padding (1 on each side)
        available_text_width = self.current_width - panel_borders_and_padding
        
        # Ensure sequence has correct length
        seq = self.reference_sequence.ljust(display_width, "?")[:display_width]
        
        # CRITICAL FIX: Only show full bases if the sequence length is EXACTLY equal to 
        # or less than the available text width - no flexibility
        use_condensed = (display_width > available_text_width)
        
        if use_condensed:
            # Calculate scale factor based on exact division
            scale_factor = math.ceil(display_width / available_text_width)
            
            # Ensure integer scaling
            scale_factor = max(1, int(scale_factor))
            
            condensed_seq = []
            
            # Create condensed representation
            for i in range(0, display_width, scale_factor):
                segment = seq[i:i+scale_factor]
                if len(segment) > 0:
                    if len(set(segment)) == 1:
                        # All same base
                        condensed_seq.append(segment[0])
                    else:
                        # Mixed bases - use complexity indicator
                        bases_set = set(segment.upper())
                        if 'N' in bases_set or '?' in bases_set:
                            condensed_seq.append('?')  # Unknown
                        elif len(bases_set) > 2:
                            condensed_seq.append('■')  # High complexity
                        else:
                            condensed_seq.append('▪')  # Low complexity
            
            # Create colored text
            colored_seq = Text()
            for base in condensed_seq:
                if base in '■▪':
                    style = get_color(self.theme, "HIGH_COMPLEXITY" if base == '■' else "LOW_COMPLEXITY")
                else:
                    style = get_color(self.theme, base.upper(), "default")
                colored_seq.append(base, style=style)
                
        else:
            # Full base display for smaller regions - ONLY when ALL bases will fit
            colored_seq = Text()
            for base in seq.upper():
                style = get_color(self.theme, base, "default")
                colored_seq.append(base, style=style)
        
        # Return as a panel with the same padding as other tracks
        title = "Reference"
        if use_condensed:
            title += f" (Condensed {scale_factor}:1 ratio)"
        else:
            title += " (Full resolution)"
        
        return Panel(colored_seq, title=title, 
                    border_style=get_color(self.theme, "SEQUENCE_BORDER"), 
                    box=box.ROUNDED)
        
    def _get_ruler_track(self) -> Text:
        """
        Create a position ruler with proper scaling.
        
        Returns:
            Text: Formatted ruler text
        """
        logging.debug("Getting ruler track.")
        display_width = self.end - self.start + 1
        term_width = self.current_width
        ruler = Text(style="dim")
        start_str = format_position(self.start)
        end_str = format_position(self.end)

        # Add start position
        ruler.append(start_str + " ")

        # Calculate space for ticks
        content_width = term_width - len(start_str) - len(end_str) - 4

        # Draw ticks if space allows
        if content_width > 10 and display_width > 0:
            bp_per_char = display_width / content_width
            ticks = ""
            for i in range(content_width):
                pos_marker_bp_offset = int((i + 0.5) * bp_per_char)
                current_pos = self.start + pos_marker_bp_offset
                if current_pos % 10 == 0 and bp_per_char < 5:
                    ticks += "|"  # Major tick
                elif current_pos % 5 == 0 and bp_per_char < 3:
                    ticks += "+"  # Minor tick
                else:
                    ticks += "-"  # Default line
            ruler.append(ticks)
        elif content_width > 0:
            ruler.append("-" * content_width)

        # Add end position
        ruler.append(" " + end_str)
        logging.debug(f"Generated ruler track Text length: {len(ruler)}")
        return ruler
        
    def _get_overview_track(self) -> Text:
        """
        Create the chromosome overview track.
        
        Returns:
            Text: Formatted overview text
        """
        logging.debug("Getting overview track.")
        if not self.fasta_file:
            return Text(f"[ FASTA file not loaded ]", 
                       style=get_color(self.theme, "WARNING_STYLE"))

        try:
            chrom_len = self.fasta_file.get_reference_length(self.chrom)
        except KeyError:
            return Text(f"[ Chromosome '{self.chrom}' not found ]", 
                       style=get_color(self.theme, "WARNING_STYLE"))

        if chrom_len <= 0:
            return Text(f"[ Chromosome '{self.chrom}' length {chrom_len} ]", 
                       style=get_color(self.theme, "WARNING_STYLE"))

        # Calculate position in chromosome
        start_pos = max(1, self.start)
        end_pos = min(chrom_len, self.end)
        if start_pos > end_pos:
            end_pos = start_pos

        # Calculate fractional positions
        relative_start = (start_pos - 1) / chrom_len
        relative_end = end_pos / chrom_len

        # Determine bar width
        bar_width = max(5, self.current_width - 4)

        # Calculate marker position
        marker_start_index = int(relative_start * bar_width)
        marker_end_index = int(relative_end * bar_width)

        # Clamp indices
        marker_start_index = max(0, min(bar_width, marker_start_index))
        marker_end_index = max(marker_start_index, min(bar_width, marker_end_index))

        # Ensure marker width is at least 1
        if marker_start_index == marker_end_index and start_pos <= end_pos:
            marker_width = 1
        else:
            marker_width = marker_end_index - marker_start_index

        if marker_width < 1 and (self.end >= self.start):
            marker_width = 1

        # Calculate prefix/suffix lengths
        marker_width = min(marker_width, bar_width - marker_start_index)
        prefix_len = marker_start_index
        suffix_len = bar_width - marker_start_index - marker_width

        # Final sanity checks
        prefix_len = max(0, prefix_len)
        suffix_len = max(0, suffix_len)
        marker_width = max(0, marker_width)

        # Construct bar
        overview_bar = "-" * prefix_len + "■" * marker_width + "-" * suffix_len
        result = Text(f"[{overview_bar}]", style=get_color(self.theme, "OVERVIEW_STYLE"))
        logging.debug(f"Generated overview track Text length: {len(result)}")
        return result
    
    def _get_variants_track(self) -> Panel:
        """
        Generate variants track with proper scaling and detail levels.
        
        Returns:
            Panel: Formatted variants panel
        """
        logging.debug("Getting variants track.")
        if not self.vcf_file:
            return Panel(Text(""), title="Variants", 
                        border_style=get_color(self.theme, "VARIANTS_BORDER"), 
                        box=box.ROUNDED)
        
        # Fetch variants for current region
        variants = self._fetch_variants()
        
        if not variants:
            return Panel(Text("No variants in region", style="dim"), 
                         title="Variants", 
                         border_style=get_color(self.theme, "VARIANTS_BORDER"), 
                         box=box.ROUNDED)
        
        # Determine display strategy based on zoom level
        display_width = self.end - self.start + 1
        term_width = self.current_width - 10  # Account for panel padding
        show_detailed = display_width <= VARIANT_DISPLAY_THRESHOLD
        show_details_list = display_width <= DETAILED_VIEW_THRESHOLD
        
        # Map variant positions to screen coordinates
        view_start = self.start - 1  # Convert to 0-based
        
        # Initialize display
        marker_display = [' '] * term_width
        marker_styles = ['default'] * term_width
        details = []
        
        # Process each variant
        for pos, variant in sorted(variants.items()):
            # Convert to screen coordinates
            screen_pos = genomic_to_screen_coord(pos - 1, view_start, display_width, term_width)
            
            # Only show if within screen bounds
            if 0 <= screen_pos < term_width:
                # Get variant type
                ref = variant.get('ref', '.')
                alt = variant.get('alt', '.').split(',')[0]  # Take first alt allele
                
                # Determine symbol and style based on variant type and zoom level
                if len(alt) == 1 and len(ref) == 1:  # SNP
                    if show_detailed:
                        # Show actual ALT base when zoomed in
                        marker_display[screen_pos] = alt
                        marker_styles[screen_pos] = get_color(self.theme, alt.upper(), "default")
                    else:
                        # Show caret for SNPs when zoomed out
                        marker_display[screen_pos] = '^'
                        marker_styles[screen_pos] = get_color(self.theme, "VARIANTS_BORDER")
                elif len(alt) > len(ref):  # Insertion
                    marker_display[screen_pos] = '+'
                    marker_styles[screen_pos] = get_color(self.theme, "INSERTION_STYLE")
                elif len(alt) < len(ref):  # Deletion
                    marker_display[screen_pos] = '-'
                    marker_styles[screen_pos] = get_color(self.theme, "DELETION_STYLE")
                else:  # Complex substitution
                    marker_display[screen_pos] = '*'
                    marker_styles[screen_pos] = get_color(self.theme, "MISMATCH_STYLE")
            
            # Always add details for list display
            detail = f"{format_position(pos)}: {variant.get('ref', '.')}>{variant.get('alt', '.')}"
            
            # Add additional details if showing details list
            if show_details_list:
                attrs = []
                if 'qual' in variant:
                    attrs.append(f"QUAL={variant['qual']:.1f}")
                for key in ['DP', 'AF', 'AC']:
                    if key in variant:
                        attrs.append(f"{key}={variant[key]}")
                if attrs:
                    detail += " [" + ", ".join(attrs) + "]"
                    
            details.append(detail)
        
        # Create display content
        content = Text()
        
        # Add marker line with proper styling for each marker
        marker_line = Text()
        for char, style in zip(marker_display, marker_styles):
            marker_line.append(char, style=style)
        content.append(marker_line)
        
        # Add details
        content.append("\n\n")
        content.append(Text("\n".join(details[:10]), style=get_color(self.theme, "VARIANTS_BORDER")))
        if len(details) > 10:
            content.append("\n... and " + str(len(details) - 10) + " more variants")
        
        # Create panel title with variant count and display mode info
        title = f"Variants ({len(variants)} found)"
        if show_detailed:
            title += " - Showing actual bases"
        
        return Panel(content, title=title, 
                    border_style=get_color(self.theme, "VARIANTS_BORDER"), 
                    box=box.ROUNDED)
    
    def _get_coverage_track(self) -> Panel:
        """
        Generate coverage track with proper scaling.
        
        Returns:
            Panel: Formatted coverage panel
        """
        if not self.bam_file:
            return Panel(Text(""), title="Coverage", 
                        border_style=get_color(self.theme, "COVERAGE_BORDER"), 
                        box=box.ROUNDED)
        
        display_width = self.end - self.start + 1
        term_width = self.current_width - 10  # Account for panel padding
        view_start = self.start - 1  # Convert to 0-based
        
        try:
            # Process alignments
            alignments = self._process_alignments()
            
            if not alignments:
                return Panel(Text("No alignments in region", style="dim"), 
                             title="Coverage", 
                             border_style=get_color(self.theme, "COVERAGE_BORDER"), 
                             box=box.ROUNDED)
            
            # Calculate coverage
            coverage = self._calculate_coverage(alignments, view_start, display_width)
            
            # Scale coverage to screen width
            scaled_coverage = self._scale_coverage(coverage, display_width, term_width)
            
            # Generate display rows
            display_rows, max_coverage = self._get_coverage_display(scaled_coverage)
            
            # Convert to Text
            content = Text()
            for row in display_rows:
                row_text = Text()
                for char, style in row:
                    row_text.append(char, style=style)
                content.append(row_text)
                content.append("\n")
            
            # Add scale indicator
            content.append(f"Max depth: {int(max_coverage)}", style="dim cyan")
            
            return Panel(content, title="Coverage", 
                        border_style=get_color(self.theme, "COVERAGE_BORDER"), 
                        box=box.ROUNDED)
            
        except Exception as e:
            logging.exception("Error generating coverage track")
            return Panel(Text(f"Error: {e}", style=get_color(self.theme, "ERROR_STYLE")), 
                         title="Coverage", 
                         border_style=get_color(self.theme, "ERROR_STYLE"), 
                         box=box.ROUNDED)
    
    def _get_alignments_track(self) -> Panel:
        """
        Generate alignments track with proper scaling and layout.
        
        Returns:
            Panel: Formatted alignments panel
        """
        logging.debug(f"Getting alignments track")
        if not self.bam_file:
            return Panel(Text(""), title="Alignments", 
                        border_style=get_color(self.theme, "ALIGNMENTS_BORDER"), 
                        box=box.ROUNDED)
        
        view_width = self.end - self.start + 1
        view_start = self.start - 1  # Convert to 0-based
        term_width = self.current_width - 10  # Account for panel padding
        
        # Skip display if region is extremely large
        if view_width > MAX_REGION_WIDTH * 2:
            content = Text(f"Region too large ({format_position(view_width)} bp) for alignment display. " +
                          f"Zoom in for detailed view.", 
                          style=get_color(self.theme, "WARNING_STYLE"))
            return Panel(content, title=f"Alignments (region too large)", 
                        border_style=get_color(self.theme, "ALIGNMENTS_BORDER"), 
                        box=box.ROUNDED)
        
        try:
            # Process alignments
            alignments = self._process_alignments()
            
            if not alignments:
                content = Text("No alignments in region", style="dim")
                return Panel(content, title="Alignments", 
                            border_style=get_color(self.theme, "ALIGNMENTS_BORDER"), 
                            box=box.ROUNDED)
            
            # Layout alignments into rows with proper scaling
            display_rows = self._layout_alignments(
                alignments, view_start, view_width, term_width, MAX_ALIGNMENT_ROWS
            )
            
            # Convert to Text for display
            content = Text()
            for row in display_rows:
                row_text = Text()
                for char, style in row:
                    row_text.append(char, style=style)
                content.append(row_text)
                content.append("\n")
            
            # Add title with tag info
            title = f"Alignments ({len(alignments)} reads, showing max {MAX_ALIGNMENT_ROWS} rows)"
            if self.tag_name:
                title += f" - Colored by '{self.tag_name}' tag"
            
            return Panel(content, title=title, 
                        border_style=get_color(self.theme, "ALIGNMENTS_BORDER"), 
                        box=box.ROUNDED)
            
        except Exception as e:
            logging.exception("Error processing alignments")
            return Panel(Text(f"Error: {e}", style=get_color(self.theme, "ERROR_STYLE")), 
                         title="Alignments", 
                         border_style=get_color(self.theme, "ERROR_STYLE"), 
                         box=box.ROUNDED)
    
    def update_layout(self):
        """Update all layout components with current state and proper scaling."""
        logging.debug("--- Updating Layout ---")
        self.current_width = self.console.width
        view_width = self.end - self.start + 1
        logging.debug(f"Terminal width: {self.current_width}, View width: {view_width}bp")

        # Update Header
        header_text = Text(f"Region: {self.chrom}:{format_position(self.start)}-{format_position(self.end)} "
                          f"(Size: {format_position(view_width)} bp)", justify="center")
        self.layout["header"].update(Panel(
            header_text, 
            style=get_color(self.theme, "HEADER_STYLE"), 
            box=box.ROUNDED)
        )

        # Update tracks
        self.layout["overview"].update(self._get_overview_track())
        self.layout["sequence"].update(self._get_sequence_track())
        self.layout["ruler"].update(self._get_ruler_track())

        # Update variants track
        variants_layout = self.layout.get("variants")
        if variants_layout:
            variants_layout.visible = bool(self.vcf_file)
            if variants_layout.visible:
                variants_layout.update(self._get_variants_track())

        # Update coverage track
        coverage_layout = self.layout.get("coverage")
        if coverage_layout:
            coverage_layout.visible = bool(self.bam_file)
            if coverage_layout.visible:
                coverage_layout.update(self._get_coverage_track())

        # Update alignments track
        alignments_layout = self.layout.get("alignments")
        if alignments_layout:
            alignments_layout.visible = bool(self.bam_file)
            if alignments_layout.visible:
                alignments_layout.update(self._get_alignments_track())

        # Update Footer
        display_mode = ""
        if view_width > SEQUENCE_DISPLAY_THRESHOLD:
            display_mode = f" | Condensed mode ({format_position(view_width)} bp)"
        elif view_width <= DETAILED_VIEW_THRESHOLD:
            display_mode = f" | Detailed mode ({format_position(view_width)} bp)"
            
        theme_mode = f"[t] = {'dark' if self.light_mode else 'light'} mode | "
        footer_text = f"Commands: [a] = go left | [d] = go right | [w] = zoom in | [s] = zoom out | [g]oto | {theme_mode}[q]uit{display_mode}"
        self.layout["footer"].update(Panel(
            Text(footer_text, justify="center"), 
            style=get_color(self.theme, "FOOTER_STYLE"), 
            box=box.ROUNDED)
        )
        logging.debug("--- Layout Update Finished ---")

    def run(self):
        """Run the main interactive loop."""
        logging.info("Starting main run loop...")
        with Live(self.layout, console=self.console, refresh_per_second=4, vertical_overflow="visible") as live:
            while True:
                # Update layout and refresh
                try:
                    self.update_layout()
                    live.refresh()
                except Exception as e:
                    logging.exception("Error during layout update")
                    live.stop()
                    self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error during layout update:[/]")
                    self.console.print_exception(show_locals=True)
                    print("\nAn unexpected error occurred. Exiting.", file=sys.stderr)
                    sys.exit(1)

                # Get user command
                try:
                    live.stop()
                    command = self.console.input(f"[{get_color(self.theme, 'PROMPT_STYLE')}]Enter command (w/a/s/d/g/t/q):[/] ").strip().lower()
                    live.start()

                except KeyboardInterrupt:
                    logging.info("KeyboardInterrupt received")
                    live.stop()
                    break
                except EOFError:
                    logging.info("EOFError received")
                    live.stop()
                    break
                except Exception as e:
                    logging.exception("Error during input")
                    live.stop()
                    self.console.print(f"[{get_color(self.theme, 'ERROR_STYLE')}]Error during input:[/]")
                    print("\nTrying to continue...", file=sys.stderr)
                    time.sleep(1)
                    try:
                        live.start()
                        continue
                    except:
                        print("\nFailed to recover. Exiting.", file=sys.stderr)
                        sys.exit(1)

                # Process command
                if command == 'q':
                    live.stop()
                    break
                    
                elif command == 'd':
                    self.navigate(forward=True)
                    
                elif command == 'a':
                    self.navigate(forward=False)
                    
                elif command == 'w':
                    self.zoom(zoom_in=True)
                    
                elif command == 's':
                    self.zoom(zoom_in=False)
                    
                elif command == 't':
                    self.toggle_theme()
                    
                elif command == 'g':
                    try:
                        live.stop()
                        new_region_str = self.console.input(f"[{get_color(self.theme, 'PROMPT_STYLE')}]Enter region (e.g., chr1:1000-2000):[/] ").strip()
                        live.start()
                        
                        parsed = parse_region(new_region_str)
                        if parsed:
                            if not self.set_region(*parsed):
                                live.stop()
                                self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Invalid region. Please try again.[/]")
                                time.sleep(1.5)
                                live.start()
                        else:
                            live.stop()
                            self.console.print(f"[{get_color(self.theme, 'WARNING_STYLE')}]Could not parse region. Format: chr:start-end[/]")
                            time.sleep(1.5)
                            live.start()
                            
                    except KeyboardInterrupt:
                        live.stop()
                        self.console.print(f"\n[{get_color(self.theme, 'WARNING_STYLE')}]Goto cancelled.[/]")
                        time.sleep(1)
                        live.start()
                    except EOFError:
                        live.stop()
                        break