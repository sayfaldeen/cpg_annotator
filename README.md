# CpG Site Annotator

This script annotates CpG sites with corresponding gene names for different array types (EPICv1, EPICv2, and MSA).

## Features

- Supports EPICv1, EPICv2, and MSA array types
- Automatically downloads annotation files if not provided
- Provides gene annotations including gene names, chromosomes, and positions
- Can save results to a file or return as a DataFrame

## Installation

1. Clone this repository
2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

Alternatively, you can run a `pip install .` in this directory with the setup.py and install into env

## Usage

The script can be used from the command line or imported as a Python module.

### Command Line Usage

```bash
python cpg-annotator.py input_file.txt array_type [--annotation_file path/to/annotation.tsv.gz] [--output_file results.tsv]
```

Arguments:
- `input_file.txt`: Text file containing one CpG site ID per line
- `array_type`: Type of array ('EPICv1', 'EPICv2', or 'MSA')
- `--annotation_file`: (Optional) Path to custom annotation file
- `--output_file`: (Optional) Path to save the annotated results

Example:
```bash
python cpg-annotator.py cpg_sites.txt EPICv1 --output_file annotated_results.tsv
```

### Python Module Usage

```python
from cpg_annotator import CpGAnnotator

# Initialize annotator
annotator = CpGAnnotator(array_type='EPICv1')

# List of CpG sites to annotate
cpg_list = ['cg00000029', 'cg00000108', 'cg00000109']

# Annotate CpG sites
results = annotator.annotate_cpg_sites(cpg_list)

# Save results to file
results.to_csv('annotated_results.tsv', sep='\t', index=False)
```

## Input File Format

The input file should be a text file with one CpG site ID per line:
```
cg00000029
cg00000108
cg00000109
```

## Output Format

The output is a tab-separated file with the following columns:
- probeID: The CpG site ID
- gene: The associated gene name
- chromosome: The chromosome location
- position: The genomic position 
