#!/usr/bin/env python3
"""
Script to annotate CpG sites with corresponding gene names.
Supports EPICv1, EPICv2, and MSA array types.
Can use local annotation files or download them if not provided.
"""

import logging
import os
import requests
from pathlib import Path
from typing import Optional
import polars as pl
from tqdm import tqdm
import argparse
from dataclasses import dataclass
from typing import Dict, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class AnnotationConfig:
    """Configuration for annotation settings."""
    chunk_size: int = 100000
    output_format: str = 'tsv'
    verbose: bool = False
    verify_downloads: bool = True

class CpGAnnotator:
    """A class to annotate CpG sites with gene information."""
    
    SUPPORTED_ARRAYS = ["EPICV1", "EPICV2", "MSA"]
    REQUIRED_COLUMNS = ["Probe_ID", "CpG_chrm", "CpG_beg", "CpG_end"]
    
    def __init__(self, array_type: str, config: Optional[AnnotationConfig] = None):
        """Initialize the CpG annotator.
        
        Args:
            array_type: Type of array ('EPICv1', 'EPICv2', or 'MSA')
            config: Optional configuration settings
        """
        self.array_type = array_type.upper()
        if self.array_type not in self.SUPPORTED_ARRAYS:
            raise ValueError(f"Unsupported array type: {array_type}. Must be one of {', '.join(self.SUPPORTED_ARRAYS)}")
        
        self.annotation_data = None
        self.config = config or AnnotationConfig()
        
        if self.config.verbose:
            logger.setLevel(logging.DEBUG)

    def get_annotation_url(self, array_type: str) -> str:
        """Get the URL for downloading the annotation file for the specified array type."""
        base_url = "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno"
        urls = {
            "EPICV1": f"{base_url}/EPIC/EPIC.hg38.manifest.tsv.gz",
            "EPICV2": f"{base_url}/EPICv2/EPICv2.hg38.manifest.tsv.gz",
            "MSA": f"{base_url}/MSA/MSA.hg38.manifest.tsv.gz"
        }
        return urls.get(array_type.upper())

    def download_annotation_file(self, output_dir: str = ".") -> str:
        """Download the annotation file for the specified array type.
        
        Args:
            output_dir: Directory to save the downloaded file
            
        Returns:
            Path to the downloaded file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        url = self.get_annotation_url(self.array_type)
        output_file = output_dir / f"{self.array_type.lower()}_annotation.tsv.gz"
        
        # Check if file already exists and is valid
        if output_file.exists() and self.config.verify_downloads:
            logger.info(f"File {output_file} already exists. Skipping download.")
            return str(output_file)
        
        logger.info(f"Downloading annotation file from {url}")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_file, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, desc="Downloading") as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
        
        return str(output_file)

    def load_annotation_data(self, annotation_file: str) -> None:
        """
        Load annotation data from a file.
        
        Args:
            annotation_file (str): Path to the annotation file
        """
        self.annotation_data = pl.read_csv(
            annotation_file,
            separator='\t',
            has_header=True,
            infer_schema_length=None
        )
        
        # Ensure required columns are present
        missing_columns = [col for col in self.REQUIRED_COLUMNS if col not in self.annotation_data.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns in annotation file: {', '.join(missing_columns)}")

    def annotate_cpg_sites(self, cpg_list: list[str], annotation_file: str | None = None, output_file: str | None = None) -> pl.DataFrame:
        """Annotate a list of CpG sites with gene information.
        
        Args:
            cpg_list: List of CpG site IDs to annotate
            annotation_file: Optional path to local annotation file
            output_file: Optional path to save results
            
        Returns:
            DataFrame with annotated CpG sites
        """
        if not cpg_list:
            raise ValueError("Empty list of CpG sites provided")
        
        # Load annotation data if needed
        if annotation_file:
            self.load_annotation_data(annotation_file)
        elif self.annotation_data is None:
            raise ValueError("No annotation data loaded. Call load_annotation_data() first or provide annotation_file.")
        
        # Create input DataFrame
        input_df = pl.DataFrame({"Probe_ID": cpg_list})
        
        # Process in chunks if needed
        if len(cpg_list) > self.config.chunk_size:
            logger.info(f"Processing {len(cpg_list)} CpG sites in chunks of {self.config.chunk_size}")
            results = []
            for i in range(0, len(cpg_list), self.config.chunk_size):
                chunk = cpg_list[i:i + self.config.chunk_size]
                chunk_df = pl.DataFrame({"Probe_ID": chunk})
                chunk_result = chunk_df.join(
                    self.annotation_data,
                    on="Probe_ID",
                    how="left"
                )
                results.append(chunk_result)
            results_df = pl.concat(results)
        else:
            results_df = input_df.join(
                self.annotation_data,
                on="Probe_ID",
                how="left"
            )
        
        # Save results if output file specified
        if output_file:
            output_path = Path(output_file)
            if self.config.output_format == 'csv':
                results_df.write_csv(output_path)
            else:
                results_df.write_csv(output_path, separator='\t')
            logger.info(f"Results saved to {output_file}")
        
        return results_df

def main():
    """Main function to demonstrate usage."""
    parser = argparse.ArgumentParser(
        description='Annotate CpG sites with gene information',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    required_group = parser.add_argument_group('Required Arguments')
    required_group.add_argument('input_file', help='Input file containing list of CpG sites')
    required_group.add_argument('array_type', choices=['EPICv1', 'EPICv2', 'MSA'], help='Type of array')
    
    # Input/Output options
    io_group = parser.add_argument_group('Input/Output Options')
    io_group.add_argument('--annotation_file', help='Path to local annotation file (optional)')
    io_group.add_argument('--output_file', help='Path to save annotated results (optional)')
    io_group.add_argument('--format', choices=['tsv', 'csv'], default='tsv', help='Output format')
    
    # Processing options
    processing_group = parser.add_argument_group('Processing Options')
    processing_group.add_argument('--chunk-size', type=int, default=100000, help='Process data in chunks of this size')
    processing_group.add_argument('--no-verify', action='store_true', help='Skip verification of downloaded files')
    
    # Logging options
    logging_group = parser.add_argument_group('Logging Options')
    logging_group.add_argument('--verbose', '-v', action='count', help='Increase verbosity level')
    logging_group.add_argument('--quiet', '-q', action='store_true', help='Suppress all output except errors')
    
    # Version information
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    
    args = parser.parse_args()
    
    # Configure logging based on verbosity and quiet options
    if args.quiet:
        log_level = logging.ERROR
    else:
        log_level = logging.INFO - (args.verbose * 10)
    logging.getLogger().setLevel(max(logging.DEBUG, log_level))
    
    # Create configuration
    config = AnnotationConfig(
        chunk_size=args.chunk_size,
        output_format=args.format,
        verbose=args.verbose > 0,
        verify_downloads=not args.no_verify
    )
    
    try:
        # Read CpG sites from input file
        with open(args.input_file, 'r') as f:
            cpg_list = [line.strip() for line in f if line.strip()]
        
        if not cpg_list:
            logger.error(f"No valid CpG sites found in {args.input_file}")
            return 1
        
        # Initialize annotator
        annotator = CpGAnnotator(args.array_type, config)
        
        # Annotate CpG sites
        results = annotator.annotate_cpg_sites(
            cpg_list,
            annotation_file=args.annotation_file,
            output_file=args.output_file
        )
        
        # Print summary
        logger.info(f"Annotated {len(results)} CpG sites")
        logger.info(f"Found {results['gene'].null_count()} sites with gene annotations")
        return 0
        
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

if __name__ == "__main__":
    # Call the main() function and use its return value (0 for success, 1 for error)
    # as the program's exit code when run as a script
    exit(main())