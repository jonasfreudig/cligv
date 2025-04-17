"""
Alignment representation module for clIGV.

This module contains the AlignmentRead class that represents and formats 
BAM read alignments for display in the terminal.
"""

from typing import List, Tuple


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


class AlignmentRead:
    """
    Class to store and display read alignment data with proper scaling.
    
    This class represents a single read alignment from a BAM file,
    and handles the display formatting with proper scaling for the terminal.
    """
    
    def __init__(self, read, ref_start, ref_end, has_tag=False, tag_value=None, tag_style=None):
        """
        Initialize an alignment read.
        
        Args:
            read: The pysam read object
            ref_start: 0-based reference start position
            ref_end: 0-based reference end position
            has_tag: Whether the read has the specified tag
            tag_value: Value of the tag (if has_tag is True)
            tag_style: Rich style string for the tag
        """
        self.read = read
        self.ref_start = ref_start  # 0-based reference start position
        self.ref_end = ref_end      # 0-based reference end position
        self.query_name = read.query_name
        self.is_reverse = read.is_reverse
        self.has_tag = has_tag
        self.tag_value = tag_value
        self.tag_style = tag_style
        self.representation = []    # List of (pos, char, style) tuples
        
    def add_position(self, ref_pos: int, char: str, style: str) -> None:
        """
        Add a position to the read representation.
        
        Args:
            ref_pos: Reference position (0-based)
            char: Character to display at this position
            style: Rich style string for this character
        """
        self.representation.append((ref_pos, char, style))
        
    def covers_position(self, pos: int) -> bool:
        """
        Check if read covers a specific position.
        
        Args:
            pos: Genomic position to check (0-based)
            
        Returns:
            bool: True if the read covers this position
        """
        return self.ref_start <= pos <= self.ref_end
        
    def generate_display(self, view_start: int, view_width: int, 
                         screen_width: int) -> List[Tuple[str, str]]:
        """
        Generate display representation for this read with proper scaling.
        
        Args:
            view_start: Start of current view (0-based)
            view_width: Width of genomic region being viewed
            screen_width: Available screen width for display
            
        Returns:
            List of (char, style) tuples for screen display
        """
        # Sort representation by position
        self.representation.sort(key=lambda x: x[0])
        
        # Initialize display array with spaces
        display = [(' ', 'default')] * screen_width
        
        # Map each genomic position to screen position
        for ref_pos, char, style in self.representation:
            screen_pos = genomic_to_screen_coord(ref_pos, view_start, view_width, screen_width)
            
            # Only update if position is within screen bounds
            if 0 <= screen_pos < screen_width:
                current_char = display[screen_pos][0]
                
                # Priority for display: more specific information overwrites less specific
                if current_char == ' ' or current_char == '.' or char not in (' ', '.'):
                    display[screen_pos] = (char, style)
        
        return display