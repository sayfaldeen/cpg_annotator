"""
Tests for the CpG annotator script using real annotation files and CpG sites.
"""

import os
import pytest
import tempfile
import shutil
from pathlib import Path
import sys

# Add the parent directory to the Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cpg_annotator.annotator import CpGAnnotator

# Test data
TEST_CPG_FILE = "tests/data/test_cpg_sites.txt"
TEST_ANNOTATION_FILE = "tests/data/test_annotation.tsv.gz"

@pytest.fixture
def test_annotator():
    """Create a test annotator instance."""
    return CpGAnnotator("EPICv1")

@pytest.fixture
def download_test_annotation(tmp_path):
    """Download test annotation file."""
    annotator = CpGAnnotator("EPICv1")
    return annotator.download_annotation_file(output_dir=str(tmp_path))

@pytest.fixture
def test_cpg_file(tmp_path):
    """Create a test CpG sites file."""
    cpg_file = tmp_path / "test_cpg_sites.txt"
    with open(cpg_file, "w") as f:
        f.write("cg00000029\ncg00000108\ncg00000109\n")
    return str(cpg_file)

def test_initialization(test_annotator):
    """Test annotator initialization."""
    assert test_annotator.array_type == "EPICV1"
    assert test_annotator.annotation_data is None

def test_get_annotation_url(test_annotator):
    """Test getting annotation URL."""
    url = test_annotator.get_annotation_url("EPICv1")
    assert "EPIC.hg38.manifest.tsv.gz" in url

def test_load_annotation_data(test_annotator, download_test_annotation):
    """Test loading annotation data."""
    test_annotator.load_annotation_data(download_test_annotation)
    assert test_annotator.annotation_data is not None
    assert "Probe_ID" in test_annotator.annotation_data.columns

def test_annotate_cpg_sites(test_annotator, download_test_annotation):
    """Test annotating CpG sites with gene information."""
    # Read test CpG sites
    with open(TEST_CPG_FILE, 'r') as f:
        cpg_list = [line.strip() for line in f if line.strip()]

    # Test annotation
    results = test_annotator.annotate_cpg_sites(
        cpg_list,
        annotation_file=download_test_annotation
    )

    assert len(results) == len(cpg_list)
    assert "Probe_ID" in results.columns
    assert "CpG_chrm" in results.columns
    assert results["Probe_ID"].to_list() == cpg_list

def test_output_file(test_annotator, download_test_annotation, tmp_path):
    """Test saving results to file."""
    output_file = tmp_path / "results.tsv"
    with open(TEST_CPG_FILE, 'r') as f:
        cpg_list = [line.strip() for line in f if line.strip()]

    results = test_annotator.annotate_cpg_sites(
        cpg_list,
        annotation_file=download_test_annotation,
        output_file=str(output_file)
    )

    assert output_file.exists()
    assert len(results) == len(cpg_list)

def test_invalid_cpg_sites(test_annotator, download_test_annotation):
    """Test handling of invalid CpG sites."""
    invalid_cpg_list = ['invalid_cpg_1', 'invalid_cpg_2']
    results = test_annotator.annotate_cpg_sites(
        invalid_cpg_list,
        annotation_file=download_test_annotation
    )

    # Verify that invalid sites are handled gracefully
    assert len(results) == len(invalid_cpg_list)
    assert results["Probe_ID"].to_list() == invalid_cpg_list
    assert results["CpG_chrm"].is_null().all()

def test_different_array_types():
    """Test different array types."""
    for array_type in ['EPICv1', 'EPICv2', 'MSA']:
        annotator = CpGAnnotator(array_type)
        assert annotator.array_type == array_type.upper()

def test_cleanup(tmp_path):
    """Test cleanup of temporary files."""
    annotator = CpGAnnotator("EPICv1")
    annotation_file = annotator.download_annotation_file(output_dir=str(tmp_path))
    
    # Verify file exists
    assert os.path.exists(annotation_file)
    
    # Clean up
    shutil.rmtree(tmp_path)
    
    # Verify file is gone
    assert not os.path.exists(annotation_file) 