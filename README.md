# Chloroplast Genome Analysis Suite

[![License: MIT] (https://opensource.org/license/mit)
[![Python 3.7+]	(https://www.python.org/downloads/)
[![Version](https://github.com/abdullah30/chloroplast-genome-scripts)

A comprehensive Python toolkit for chloroplast genome analysis with 7 specialized modes. Analyze GenBank files and FASTA alignments to generate publication-ready results.

**Choose your preferred method:** Unified Script | Python Module | Individual Scripts | Jupyter Notebook | Linux CLI

---

## Table of Contents

- [Features](#-features)
- [Quick Start](#-quick-start)
- [Installation](#-installation)
- [Usage Methods](#-usage-methods)
- [Analysis Modes](#-analysis-modes)
- [Examples](#-examples)
- [Troubleshooting](#-troubleshooting)
- [Documentation](#-documentation)
- [Citation](#-citation)
- [License](#-license)

---

## Features

? **7 Analysis Modes** - Gene content, codon usage, amino acids, SNPs, introns, comprehensive data about genome regions and their GC content, along with GC content of different genes, gene tables  
? **Flexible Usage** - Unified script, module imports, individual scripts, Jupyter notebooks  
? **Publication-Ready** - Formatted Excel and Word outputs  
? **Smart Handling** - Correctly handles trans-spliced genes (rps12)  
? **Cross-Platform** - Windows, Linux, macOS  
? **Well-Documented** - Comprehensive guides included

---

## Quick Start

```bash
# 1. Install dependencies
pip install biopython pandas openpyxl numpy python-docx

# 2. Place your GenBank files (.gb) in the working directory

# 3. Run all analyses
python chloroplast_unified_analysis.py
```

**Done!** Check for Excel and Word files in your directory.

---

## Installation

### Windows
```cmd
pip install biopython pandas openpyxl numpy python-docx
```

### Linux/Mac
```bash
pip3 install biopython pandas openpyxl numpy python-docx
```

### Anaconda
```bash
conda install -c conda-forge biopython pandas openpyxl numpy
pip install python-docx
```

---

## Usage Methods

### Method 1: Unified Script (Easiest)

**Best for:** Running all analyses at once

```bash
python chloroplast_unified_analysis.py
```

**What it does:**
- Auto-detects GenBank and FASTA files
- Runs all applicable analyses
- Generates timestamped output files
- Shows progress and summary

---

### Method 2: Python Module

**Best for:** Programmatic control, custom workflows

```python
from chloroplast_analyzer import (
    run_mode_1,  # Gene Content
    run_mode_2,  # Gene Table (Word)
    run_mode_3,  # IR Boundaries
    run_mode_4,  # Codon Usage
    run_mode_5,  # Amino Acids
    run_mode_6,  # SNP Analysis
    run_mode_7,  # Introns
    run_all_modes
)

# Run specific mode
run_mode_1()	#Gene Content

# Run all modes
summary = run_all_modes()
```

**Command line:**
```bash
# Run all
python chloroplast_analyzer.py --all

# Run specific modes
python chloroplast_analyzer.py --mode 1 4 7

# List modes
python chloroplast_analyzer.py --list
```

---

### Method 3: Individual Scripts

**Best for:** Single analysis, testing

```bash
python "1_Genome_gene_count_and_summary.py"
python "2_Gene table final.py"
python "5_amino acid.py"
```

---

### Method 4: Jupyter Notebook

**Best for:** Interactive analysis, visualization

```python
# Cell 1: Install (run once)
!pip install biopython pandas openpyxl numpy python-docx -q

# Cell 2: Import
from chloroplast_analyzer import run_mode_1, run_all_modes

# Cell 3: Run
run_mode_1()

# Cell 4: View results
import pandas as pd
df = pd.read_excel('Chloroplast_Gene_Analysis_*.xlsx')
display(df)

```

**Or use magic command:**
```python
%run "1_Genome_gene_count_and_summary.py"
%run "7_intron analysis"	#place the require script or all scripts in working directory	
%run "6_amino acid"	#Use this pattern for each script
```

---

### Method 5: Linux Command Line

**Best for:** Server environments, automation

```bash
# Basic usage
cd ~/chloroplast_analysis
python3 chloroplast_unified_analysis.py

# Virtual environment
python3 -m venv chloroplast_env
source chloroplast_env/bin/activate
pip install biopython pandas openpyxl numpy python-docx
python chloroplast_unified_analysis.py

# Background execution
nohup python3 chloroplast_unified_analysis.py > analysis.log 2>&1 &

# Batch processing
for dir in dataset*/; do
    cd "$dir"
    python3 ../chloroplast_unified_analysis.py
    cd ..
done
```

---

## ?? Analysis Modes

| Mode | Name | Input | Output | Description |
|------|------|-------|--------|-------------|
| **1** | Gene Content | *.gb | Excel | Counts genes, identifies duplications, identifies unique genes, gives gene count for each species |
| **2** | Gene Tables | *.gb | Word | Publication-ready tables which show each category of gene |
| **3** | IR Boundaries | *.gb | Excel | Inverted repeat, LSC, and SSC analysis, GC content of all regions, GC content of functional genes |
| **4** | Codon Usage | *.gb | Excel | Codon frequency calculation in term of Relative synonymous codon usage |
| **5** | Amino Acids | *.gb | Excel | Amino acid composition |
| **6** | SNP Analysis | *.fasta | Excel | SNP detection, Ts/Tv ratios |
| **7** | Intron Analysis | *.gb | Excel | Intron positions/lengths |

### Input Requirements

**GenBank files** (Modes 1-5, 7):
- Format: `.gb`, `.gbk`, `.genbank`, `.gbff`
- From NCBI or similar sources
- Properly annotated

**FASTA alignments** (Mode 6):
- Format: `.fasta`, `.fa`
- Must be aligned (same length)
- Use MAFFT, MUSCLE, or Clustal

### Output Files

All outputs are timestamped:
- Excel: `Analysis_Name_YYYYMMDD_HHMMSS.xlsx`
- Word: `Table_[genome_name].docx`

---

## Examples

### Example 1: Complete Analysis

```bash
python chloroplast_unified_analysis.py
```

### Example 2: Selective Modes

```python
from chloroplast_analyzer import run_mode_1, run_mode_4

run_mode_1()  # Genes
run_mode_4()  # Codons
```

### Example 3: Jupyter Workflow

```python
# Setup
!pip install biopython pandas openpyxl numpy python-docx -q

# Run
from chloroplast_analyzer import run_all_modes
summary = run_all_modes()

# Visualize
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_excel('Chloroplast_Gene_Analysis_*.xlsx')
plt.bar(df['Species'], df['Total_Genes'])
plt.show()
```

### Example 4: Linux Automation

```bash
#!/bin/bash
# Batch process all directories

for species in */; do
    echo "Processing $species"
    cd "$species"
    python3 ../chloroplast_unified_analysis.py
    cd ..
done
```

---

## Troubleshooting

### Common Issues

**"No module named 'Bio'"**
```bash
pip install biopython
```

**"No module named 'docx'"**
```bash
pip install python-docx
```

**"No GenBank files found"**
- Check files are `.gb`, `.gbk`, `.genbank`, or `.gbff`
- Ensure files are in the working directory

**"Permission denied" (Linux)**
```bash
chmod +x chloroplast_unified_analysis.py
# or
python3 chloroplast_unified_analysis.py
```

## ?? Documentation

### Complete Guides

- **Jupyter usage, widgets, visualization** - 
- **LINUX_GUIDE** - Linux installation, virtual env, bash scripts
- **MODULE_GUIDE**- Module API, programmatic usage
- **INDIVIDUAL_SCRIPTS_README** - Individual script details
- **METHODS** - Academic methodology
- **LICENSE_RECOMMENDATION** - License guide

---

## Citation

If you use this software in your research, please cite:

```
Abdullah. (2025). Chloroplast Genome Analysis Suite (Version 1.0). 
GitHub: https://github.com/abdullah30/chloroplast-genome-scripts
```

**BibTeX:**
```bibtex
@software{abdullah2025chloroplast,
  author = {Abdullah},
  title = {Chloroplast Genome Analysis Suite},
  year = {2025},
  version = {1.0},
  url = {https://github.com/abdullah30/chloroplast-genome-scripts}
}
```

---

## Acknowledgments

Built with:
- **Biopython** - Sequence analysis
- **Pandas** - Data manipulation  
- **OpenPyXL** - Excel files
- **NumPy** - Numerical operations
- **python-docx** - Word documents

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

**You can:** Use, modify, distribute freely (even commercially)  
**You must:** Include license and copyright notice

---

##  Quick Command Reference

```bash
# Install
pip install biopython pandas openpyxl numpy python-docx

# Run all analyses
python chloroplast_unified_analysis.py

# Run specific mode
python chloroplast_analyzer.py --mode 1

# Jupyter
%run chloroplast_unified_analysis.py

# Linux background
nohup python3 chloroplast_unified_analysis.py > log.txt 2>&1 &

# Check results
ls -lh *.xlsx *.docx
```

---

## Requirements

**System:** Python 3.7+  
**Packages:** biopython, pandas, openpyxl, numpy, python-docx  
**Input:** GenBank files (.gb) and/or FASTA alignments (.fasta)  
**Output:** Excel (.xlsx) and Word (.docx) files

Install all:
```bash
pip install -r requirements.txt
```

---

## Contact

**Author:** Abdullah  
**GitHub:** [abdullah30/chloroplast-genome-scripts](https://github.com/abdullah30/chloroplast-genome-scripts)

---

## Contributing

Contributions welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

---

## Star This Repository

If you find this useful, please star the repository!

---

**Made with for the chloroplast genomics community**

**Version 1.0** | **MIT License** | **Python 3.7+**

---

*Last updated: December 2025*
