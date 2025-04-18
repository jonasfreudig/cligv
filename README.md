# clIGV - command line Interactive Genome Viewer

[![PyPI version](https://img.shields.io/pypi/v/cligv.svg)](https://pypi.org/project/cligv/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/cligv.svg)](https://anaconda.org/conda-forge/cligv)
[![License](https://img.shields.io/github/license/jonasfreudig/cligv.svg)](https://github.com/jonasfreudig/cligv/blob/main/LICENSE)

A fast, interactive genome browser for the command line. clIGV (command line Interactive Genome Viewer) allows you to view genomic sequences, alignments, and variants directly in your terminal with a simple interface.

<img src="docs/images/clIGV_img1.png" width="900">

## Features

- View reference genome sequences from FASTA files
- Display variant calls from VCF files
- Visualize read alignments from BAM files with coverage tracks
- Intuitive keyboard navigation and zooming
- Color-coded display with both dark and light theme options
- Tag-based read coloring for highlighting specific attributes
- Scales well for both small regions and larger genomic intervals

## Installation

### Using pip

```bash
pip install cligv
```

### Using conda

```bash
conda install -c conda-forge cligv
```

### From source

```bash
git clone https://github.com/jonasfreudig/cligv.git
cd cligv
pip install -e .
```

## Quick Start

```bash
# Basic usage with just a reference genome
cligv reference.fasta

# View with variants
cligv reference.fasta -v variants.vcf.gz

# View with read alignments
cligv reference.fasta -b alignments.bam

# Load directly to a specific region
cligv reference.fasta -r chr1:1000-2000

# Full example with all options
cligv reference.fasta -v variants.vcf.gz -b alignments.bam -r chr1:1000-2000 -t ha --light-mode
```

## Navigation

- `a` - Move left
- `d` - Move right
- `w` - Zoom in
- `s` - Zoom out
- `g` - Go to a specific region
- `t` - Toggle between dark and light mode
- `q` - Quit

## Requirements

- Python 3.8+
- pysam
- rich


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.

## Acknowledgments

- Inspired by tools like IGV, samtools tview, ASCIIGenome, and other genomic viewers
- Thanks to all contributors and users for their feedback and suggestions