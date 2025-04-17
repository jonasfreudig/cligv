"""Basic tests for cligv package."""

import pytest
from cligv.utils import parse_region, format_position


def test_format_position():
    """Test the position formatting function."""
    assert format_position(1000) == "1,000"
    assert format_position(1234567) == "1,234,567"
    assert format_position(0) == "0"


def test_parse_region():
    """Test the region parsing function."""
    # Test chromosome:start-end format
    result = parse_region("chr1:1000-2000")
    assert result == ("chr1", 1000, 2000)
    
    # Test chromosome:position format
    result = parse_region("chr2:5000")
    assert result[0] == "chr2"
    assert result[1] <= 5000
    assert result[2] >= 5000
    
    # Test chromosome-only format
    result = parse_region("chrM")
    assert result[0] == "chrM"
    assert result[1] > 0
    
    # Test with commas in numbers
    result = parse_region("chr1:1,000-2,000")
    assert result == ("chr1", 1000, 2000)
    
    # Test invalid input
    assert parse_region("invalid input") is None
    assert parse_region("") is None
    assert parse_region(None) is None