#!/usr/bin/env python3
"""Setup script for cligv package."""

from setuptools import setup

# Use information from pyproject.toml
setup(
    name="cligv",
    packages=["cligv"],
    entry_points={
        "console_scripts": [
            "cligv=cligv.cli:main",
        ],
    },
)