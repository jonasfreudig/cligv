"""
Utility functions for the command line Interactive Genome Viewer (clIGV).

This module provides various utility functions for coordinate transformations,
region parsing, and formatting.
"""

import re
import logging
from typing import Optional, Tuple


# --- Constants ---
DEFAULT_REGION_WIDTH = 250  # Initial view width
MAX_REGION_WIDTH = 100000  # Maximum width before warning
MAX_FETCH_SIZE = 1000000  # Maximum BAM fetch size


def genomic_to_screen_coord(genomic_pos: int, view_start: int, view_width: int, 
                           screen_width: int) -> int:
    """
    Convert genomic coordinate to screen coordinate with proper scaling.
    
    Args:
        genomic_pos: Genomic position (0-based)
        view_start: Start of current view (0-based)
        view_width: Width of genomic region being viewed
        screen_width: Available screen width for display
        
    Returns:
        Screen position (column) for the genomic position
    """
    # Calculate position relative to view start
    rel_pos = genomic_pos - view_start
    
    if rel_pos < 0 or rel_pos >= view_width:
        return -1  # Outside of visible range
    
    # Calculate scaling factor if region is wider than screen
    if view_width > screen_width:
        scaling_factor = view_width / screen_width
        screen_pos = int(rel_pos / scaling_factor)
    else:
        screen_pos = rel_pos
        
    return screen_pos


def screen_to_genomic_coord(screen_pos: int, view_start: int, view_width: int, 
                           screen_width: int) -> int:
    """
    Convert screen coordinate to genomic coordinate with proper scaling.
    
    Args:
        screen_pos: Screen position (column)
        view_start: Start of current view (0-based)
        view_width: Width of genomic region being viewed
        screen_width: Available screen width for display
        
    Returns:
        Genomic position (0-based)
    """
    # Calculate scaling factor if region is wider than screen
    if view_width > screen_width:
        scaling_factor = view_width / screen_width
        genomic_offset = int(screen_pos * scaling_factor)
    else:
        genomic_offset = screen_pos
        
    return view_start + genomic_offset


def format_position(pos: int) -> str:
    """Format position with commas for thousands.
    
    Args:
        pos: Position to format
        
    Returns:
        Formatted position string with thousands separators
    """
    return f"{pos:,}"


def parse_region(region_str: Optional[str]) -> Optional[Tuple[str, Optional[int], Optional[int]]]:
    """
    Parse a genomic region string (chr:start-end or chr:pos or chr).
    
    Args:
        region_str: Region string to parse
        
    Returns:
        Tuple of (chromosome, start, end) if valid, None otherwise.
        Start and end are 1-based coordinates.
    """
    if not region_str: 
        return None
        
    region_str = region_str.replace(',', '').strip()
    
    # Match Chr:Start-End
    match = re.match(r"([\w\.-]+):(\d+)-(\d+)", region_str)
    if match:
        chrom, start, end = match.groups()
        return chrom, int(start), int(end)
        
    # Match Chr:Start
    match = re.match(r"([\w\.-]+):(\d+)", region_str)
    if match:
        chrom, start = match.groups()
        start_int = int(start)
        view_start = max(1, start_int - DEFAULT_REGION_WIDTH // 2)
        view_end = view_start + DEFAULT_REGION_WIDTH - 1
        return chrom, view_start, view_end
        
    # Match Chr only
    match = re.match(r"([\w\.-]+)", region_str)
    if match:
        chrom = match.group(1)
        return chrom, 1, DEFAULT_REGION_WIDTH
        
    # Log warning if parsing fails
    logging.warning(f"Could not parse region string '{region_str}'")
    return None