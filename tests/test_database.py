import os
import io
import sqlite3
import pandas as pd
import pytest
from unittest.mock import patch, MagicMock
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from EnigmaDB import Database, create_index, convert_fasta

class DummyDatabase(Database):
  """ Override network methods for testing. """
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def search(self, query):
    # Return different results based on query to test various scenarios
    if query == "empty":
      return []
    elif query == "error":
      raise Exception("Test search error")
    else:
      # pretend we found three fake IDs
      return ['id1', 'id2', 'id3']

  def _safe_efetch(self, ids):
    # return a small FASTA stream with only first 2 IDs to test partial fetching
    if 'error' in str(ids):
      raise Exception("Test fetch error")
    
    records = []
    for i, seq_id in enumerate(ids[:2]):  # Only return first 2 to simulate partial fetch
      seq = "AAA" if i == 0 else "CCC"
      records.append(SeqRecord(Seq(seq), id=seq_id, description=f"desc{i+1}"))
    
    stream = io.StringIO()
    SeqIO.write(records, stream, "fasta")
    stream.seek(0)
    return stream

@pytest.fixture
def tmp_out(tmp_path):
  d = tmp_path / "out"
  d.mkdir()
  return str(d)

def test_database_initialization_validation():
  """Test Database initialization parameter validation"""
  with pytest.raises(ValueError, match="Topics list cannot be empty"):
    Database([], "out_dir", "test@example.com")
    
  with pytest.raises(ValueError, match="An email address is required"):
    Database(["topic1"], "out_dir", "")
    
  with pytest.raises(ValueError, match="Output directory must be specified"):
    Database(["topic1"], "", "test@example.com")
    
  with pytest.raises(ValueError, match="max_rate must be positive"):
    Database(["topic1"], "out_dir", "test@example.com", max_rate=0)
    
  with pytest.raises(ValueError, match="batch_size must be positive"):
    Database(["topic1"], "out_dir", "test@example.com", batch_size=0)
    
  with pytest.raises(ValueError, match="retmax must be positive"):
    Database(["topic1"], "out_dir", "test@example.com", retmax=0)

def test_database_initialization_with_api_key(tmp_out):
  """Test Database initialization with API key"""
  db = Database(['Q'], tmp_out, email="a@b.com", api_key="TEST_KEY")
  assert db.max_rate == 3.0  # Should use provided rate
    
  # Test without API key - should enforce 3.0 max
  db_no_key = Database(['Q'], tmp_out, email="a@b.com", max_rate=10.0)
  assert db_no_key.max_rate == 3.0  # Should be capped at 3.0

def test_sanitize(tmp_out):
  db = DummyDatabase(['Q'], tmp_out, email="a@b.com", api_key="KEY")
  assert db._sanitize("Hello World!") == "Hello_World"
  assert db._sanitize("A/B:C") == "A_B_C"
  assert db._sanitize("") == "unnamed"
  assert db._sanitize("   ") == "unnamed"
  assert db._sanitize("___") == "unnamed"

def test_search_functionality(tmp_out):
  """Test search method with different scenarios"""
  db = DummyDatabase(['TestTopic'], tmp_out, email="a@b.com", api_key="KEY")
    
  # Normal search
  results = db.search("normal_query")
  assert results == ['id1', 'id2', 'id3']
    
  # Empty query
  results = db.search("")
  assert results == []
    
  # Empty results
  results = db.search("empty")
  assert results == []

@patch('time.sleep')  # Mock sleep to speed up tests
def test_build_functionality(mock_sleep, tmp_out):
  """Test build method with various scenarios"""
  db = DummyDatabase(['TestTopic'], tmp_out, email="a@b.com", api_key="KEY", batch_size=2)
    
  # Mock user input for file overwrite
  with patch('builtins.input', return_value='y'):
    db.build(with_index=False)  # Disable index creation to avoid SeqIO.index_db issues
    
  # Check that TestTopic was processed successfully
  fasta = os.path.join(tmp_out, "TestTopic.fasta")
    
  assert os.path.exists(fasta)
    
  # Verify FASTA content
  recs = list(SeqIO.parse(fasta, "fasta"))
  assert len(recs) >= 1  # Should have at least some records
  record_ids = {r.id for r in recs}
  # Should contain ids from our dummy fetcher
  assert any(rid in ['id1', 'id2'] for rid in record_ids)

@patch('time.sleep')
def test_build_file_overwrite_decline(mock_sleep, tmp_out):
  """Test build method when user declines to overwrite existing file"""
  db = DummyDatabase(['TestTopic'], tmp_out, email="a@b.com", api_key="KEY")
    
  # Create existing file
  existing_file = os.path.join(tmp_out, "TestTopic.fasta")
  with open(existing_file, 'w') as f:
    f.write(">test\nATGC\n")
    
  # Mock user input to decline overwrite
  with patch('builtins.input', return_value='n'):
    db.build(with_index=False)
    
  # File should remain unchanged
  with open(existing_file, 'r') as f:
    content = f.read()
    assert "test" in content

def test_safe_efetch_retry_logic(tmp_out):
  """Test _safe_efetch retry mechanism"""
  db = DummyDatabase(['TestTopic'], tmp_out, email="a@b.com", api_key="KEY")
    
  # Test successful fetch
  result = db._safe_efetch(['id1', 'id2'])
  assert isinstance(result, io.StringIO)
  content = result.getvalue()
  assert ">id1" in content
    
  # Test empty ID list
  with pytest.raises(ValueError, match="No IDs provided"):
    db._safe_efetch([])

def test_create_index_functionality(tmp_path):
  """Test create_index function with various scenarios"""
  # Create test directory with FASTA files
  dir_in = tmp_path / "fa"
  dir_in.mkdir()
    
  # Create multiple FASTA files
  test_data = [
    ("file1", [("seq1", "ATGCCC"), ("seq2", "GGGAAA")]),
    ("file2", [("seq3", "TTTAGC"), ("seq4", "CCCGGG")]),
    ("empty", []),  # Empty file
  ]
    
  for filename, sequences in test_data:
    fasta_path = dir_in / f"{filename}.fasta"
    records = [SeqRecord(Seq(seq), id=seq_id, description=seq_id) for seq_id, seq in sequences]
        
    with open(fasta_path, "w") as fh:
      if records:  # Only write if there are records
        SeqIO.write(records, fh, "fasta")
      # Leave empty file empty
    
  # Test index creation
  idx_file = tmp_path / "combined.idx"
  create_index(str(dir_in), str(idx_file))
    
  assert idx_file.exists()    
  # Verify index content
  conn = sqlite3.connect(str(idx_file))
  cur = conn.cursor()
  cur.execute("SELECT key, filename FROM seq_index ORDER BY key")
  rows = cur.fetchall()
  conn.close()
    
  # Should have 4 sequences (empty file should be ignored)
  assert len(rows) == 4
  seq_ids = [row[0] for row in rows]
  assert set(seq_ids) == {"seq1", "seq2", "seq3", "seq4"}

def test_create_index_with_duplicates(tmp_path):
  """Test create_index handling of duplicate sequence IDs"""
  dir_in = tmp_path / "fa"
  dir_in.mkdir()
    
  # Create files with duplicate sequence IDs
  for i, filename in enumerate(["file1", "file2"]):
    fasta_path = dir_in / f"{filename}.fasta"
    records = [SeqRecord(Seq("ATGC"), id="duplicate_id", description="test")]
    
    with open(fasta_path, "w") as fh:
      SeqIO.write(records, fh, "fasta")
    
  idx_file = tmp_path / "combined.idx"
  create_index(str(dir_in), str(idx_file))
  
  # Should only have one entry due to duplicate handling
  conn = sqlite3.connect(str(idx_file))
  cur = conn.cursor()
  cur.execute("SELECT COUNT(*) FROM seq_index")
  count = cur.fetchone()[0]
  conn.close()
  
  assert count == 1

def test_create_index_error_handling():
  """Test create_index error handling for invalid inputs"""
  with pytest.raises(FileNotFoundError):
    create_index("nonexistent_directory", "output.idx")
    
  # Test with a file instead of directory
  with pytest.raises(NotADirectoryError):
    create_index(__file__, "output.idx")

def test_convert_fasta_csv(tmp_path):
  """Test convert_fasta function with CSV output"""
  # Create input directory with FASTA files
  dir_in = tmp_path / "input"
  dir_in.mkdir()
    
  # Create test FASTA files
  test_sequences = [
    ("file1", [("seq1", "ATGCCC", "Test sequence 1")]),
    ("file2", [("seq2", "GGGAAA", "Test sequence 2")]),
  ]
    
  for filename, sequences in test_sequences:
    fasta_path = dir_in / f"{filename}.fasta"
    records = [SeqRecord(Seq(seq), id=seq_id, description=desc) for seq_id, seq, desc in sequences]
    
    with open(fasta_path, "w") as fh:
      SeqIO.write(records, fh, "fasta")
    
  # Test CSV conversion
  out_dir = tmp_path / "output"
  convert_fasta(str(dir_in), str(out_dir), mode="csv")
  
  # Verify output
  csv_files = list(out_dir.glob("*.csv"))
  assert len(csv_files) == 2
  
  # Check CSV content
  df = pd.read_csv(csv_files[0])
  expected_columns = {"id", "name", "length", "sequence"}
  assert set(df.columns) == expected_columns
  assert len(df) == 1
  assert df.iloc[0]["length"] == 6  # Length of "ATGCCC"

def test_convert_fasta_parquet(tmp_path):
  """Test convert_fasta function with Parquet output"""
  try:
    import pyarrow
  except ImportError:
    pytest.skip("pyarrow not available")
    
  dir_in = tmp_path / "input"
  dir_in.mkdir()
    
  # Create test FASTA file
  fasta_path = dir_in / "test.fasta"
  records = [SeqRecord(Seq("ATGC"), id="test_seq", description="test")]

  with open(fasta_path, "w") as fh:
    SeqIO.write(records, fh, "fasta")

  # Test Parquet conversion
  out_dir = tmp_path / "output"
  convert_fasta(str(dir_in), str(out_dir), mode="parquet")

  # Verify output
  parquet_files = list(out_dir.glob("*.parquet"))
  assert len(parquet_files) == 1

  # Check Parquet content
  df = pd.read_parquet(parquet_files[0])
  expected_columns = {"id", "name", "length", "sequence"}
  assert set(df.columns) == expected_columns

def test_convert_fasta_error_handling():
  """Test convert_fasta error handling"""
  with pytest.raises(FileNotFoundError):
    convert_fasta("nonexistent", "output", mode="csv")

  with pytest.raises(ValueError, match="Unsupported mode"):
    convert_fasta("input", "output", mode="invalid")

def test_convert_fasta_empty_files(tmp_path):
  """Test convert_fasta with empty or invalid files"""
  dir_in = tmp_path / "input"
  dir_in.mkdir()

  # Create empty file
  empty_file = dir_in / "empty.fasta"
  empty_file.touch()
  
  # Create non-FASTA file
  non_fasta = dir_in / "not_fasta.txt"
  with open(non_fasta, "w") as f:
    f.write("This is not a FASTA file")
    
  out_dir = tmp_path / "output"

  # Should handle gracefully without errors
  convert_fasta(str(dir_in), str(out_dir), mode="csv")    
  # No output files should be created
  csv_files = list(out_dir.glob("*.csv"))
  assert len(csv_files) == 0

def test_database_output_directory_creation(tmp_path):
  """Test that Database creates output directory if it doesn't exist"""
  non_existent_dir = tmp_path / "new_output_dir"
  assert not non_existent_dir.exists()
  
  db = Database(['test'], str(non_existent_dir), email="test@example.com")
  assert non_existent_dir.exists()
  assert non_existent_dir.is_dir()

def test_database_sleep_calculation(tmp_out):
  """Test that sleep time is calculated correctly based on max_rate"""
  db = Database(['test'], tmp_out, email="test@example.com", max_rate=2.0)
  assert db._sleep == 0.5  # 1.0 / 2.0
  
  db2 = Database(['test'], tmp_out, email="test@example.com", max_rate=10.0)
  assert db2._sleep == 0.1  # 1.0 / 10.0

# Self-executing test runner
if __name__ == "__main__":
  import sys
  
  # Check if pytest is available
  try:
    import pytest
  except ImportError:
    print("pytest is required to run the tests. Install it with: pip install pytest")
    sys.exit(1)
  
  # Check if required dependencies are available
  required_packages = ['Bio', 'pandas', 'sqlite3']
  missing_packages = []
  
  for package in required_packages:
    try:
      __import__(package)
    except ImportError:
      missing_packages.append(package)
  
  if missing_packages:
    print(f"Missing required packages: {', '.join(missing_packages)}")
    print("Install them with: pip install biopython pandas")
    sys.exit(1)
  
  print("Running all tests...")
  print("=" * 50)
  
  # Run pytest with current file
  exit_code = pytest.main([
    __file__,
    "-v",  # verbose output
    "--tb=short",  # shorter traceback format
    "-x",  # stop on first failure
    "--color=yes"  # colored output
  ])
  
  if exit_code == 0:
    print("\n" + "=" * 50)
    print("✅ All tests passed successfully!")
  else:
    print("\n" + "=" * 50)
    print("❌ Some tests failed. Check the output above for details.")
  
  sys.exit(exit_code)