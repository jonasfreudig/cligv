"""
clIGV - command line Interactive Genome Viewer

A fast, interactive genome browser for the terminal that allows
visualization of genomic sequences, alignments, and variants.
"""

__version__ = "0.1.0"
__author__ = "Jonas Freudigmann"
__license__ = "MIT"

from .browser import GenomeBrowser
from .alignment import AlignmentRead