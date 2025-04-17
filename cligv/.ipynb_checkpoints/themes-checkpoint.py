"""
Color themes for the command line Interactive Genome Viewer (clIGV).

This module defines color palettes for both dark and light terminal modes,
allowing the viewer to adapt to the user's terminal preferences.
"""

# Dark theme (default) - optimized for dark backgrounds
DARK_THEME = {
    # Nucleotide colors
    "A": "bold chartreuse2",
    "C": "bold turquoise2",
    "G": "bold gold1",
    "T": "bold red",
    "N": "bold white",
    "-": "dim",
    "?": "dim",
    
    # Alignment display
    "MATCH_CHAR": ".",
    "MATCH_STYLE_FWD": "cyan",
    "MATCH_STYLE_REV": "magenta",
    "MISMATCH_STYLE": "bold white on red",
    "DELETION_CHAR": "*",
    "DELETION_STYLE": "bold magenta",
    "INSERTION_CHAR": "+",
    "INSERTION_STYLE": "bold green",
    
    # Variant display
    "VARIANT_ALLELE_STYLE": "bold white on green",
    
    # UI elements
    "HEADER_STYLE": "bold white on dark_green",
    "FOOTER_STYLE": "white on dark_green",
    "OVERVIEW_STYLE": "yellow",
    "SEQUENCE_BORDER": "dim",
    "VARIANTS_BORDER": "gold3",
    "COVERAGE_BORDER": "dark_cyan",
    "ALIGNMENTS_BORDER": "dark_cyan",
    "COVERAGE_BAR": "dark_cyan",
    "WARNING_STYLE": "yellow",
    "ERROR_STYLE": "bold red",
    "SUCCESS_STYLE": "green",
    "PROMPT_STYLE": "bold cyan",
    
    # Complexity indicators
    "HIGH_COMPLEXITY": "yellow",
    "LOW_COMPLEXITY": "cyan"
}

# Light theme - optimized for light backgrounds
LIGHT_THEME = {
    # Nucleotide colors
    "A": "bold dark_green",
    "C": "bold blue",
    "G": "bold orange3",
    "T": "bold red",
    "N": "bold black",
    "-": "dim",
    "?": "dim",
    
    # Alignment display
    "MATCH_CHAR": ".",
    "MATCH_STYLE_FWD": "dark_blue",
    "MATCH_STYLE_REV": "purple",
    "MISMATCH_STYLE": "bold white on red",
    "DELETION_CHAR": "*",
    "DELETION_STYLE": "bold magenta",
    "INSERTION_CHAR": "+",
    "INSERTION_STYLE": "bold dark_green",
    
    # Variant display
    "VARIANT_ALLELE_STYLE": "bold black on yellow",
    
    # UI elements
    "HEADER_STYLE": "bold black on green",
    "FOOTER_STYLE": "black on green",
    "OVERVIEW_STYLE": "dark_orange",
    "SEQUENCE_BORDER": "grey50",
    "VARIANTS_BORDER": "orange3",
    "COVERAGE_BORDER": "blue",
    "ALIGNMENTS_BORDER": "blue",
    "COVERAGE_BAR": "blue",
    "WARNING_STYLE": "orange3",
    "ERROR_STYLE": "bold red",
    "SUCCESS_STYLE": "dark_green",
    "PROMPT_STYLE": "bold blue",
    
    # Complexity indicators
    "HIGH_COMPLEXITY": "dark_orange",
    "LOW_COMPLEXITY": "blue"
}

def get_theme(light_mode=False):
    """Return the appropriate color theme based on mode.
    
    Args:
        light_mode (bool): If True, returns light theme for light backgrounds.
                          If False (default), returns dark theme.
    
    Returns:
        dict: Dictionary mapping theme elements to Rich-compatible color strings.
    """
    return LIGHT_THEME if light_mode else DARK_THEME

def get_color(theme, key, default=None):
    """Get a color from the theme with fallback.
    
    Args:
        theme (dict): The current theme dictionary.
        key (str): The theme color key to look up.
        default (str, optional): Default value if key not found.
    
    Returns:
        str: Rich-compatible color string.
    """
    return theme.get(key, default or "default")