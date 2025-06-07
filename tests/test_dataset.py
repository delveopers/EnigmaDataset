import torch
import pytest
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add the parent directory to the path to import EnigmaDB
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class DummyTokenizer:
  def __init__(self, mode="dna", kmer=1, continuous=False):
    self.mode = mode
    self.kmer = kmer
    self.continuous = continuous
    self.vocab_size = 4 if mode == "dna" else 20  # DNA: 4 bases, Protein: 20 amino acids
    
  def encode(self, s):
    # map A,C,G,T → 0,1,2,3 for DNA
    if self.mode == "dna":
      m = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
    else:  # protein mode
      # Simple mapping for amino acids
      amino_acids = "ACDEFGHIKLMNPQRSTVWY"
      m = {aa: i for i, aa in enumerate(amino_acids)}
    return [m.get(ch.upper(), 0) for ch in s]

# Set up mock before importing Dataset
def setup_mock():
  """Setup mock biosaic module"""
  biosaic_mock = type('biosaic', (), {})()
  biosaic_mock.Tokenizer = DummyTokenizer
  sys.modules['biosaic'] = biosaic_mock

# Apply mock immediately
setup_mock()

# Now import Dataset after mocking
from EnigmaDB import Dataset

@pytest.fixture(autouse=True)
def patch_tokenizer(monkeypatch):
  # monkey-patch biosaic.Tokenizer (note: capital T)
  import builtins
  import sys

  # Create mock biosaic module
  biosaic_mock = type('biosaic', (), {})()
  biosaic_mock.Tokenizer = DummyTokenizer

  monkeypatch.setitem(sys.modules, 'biosaic', biosaic_mock)
  yield

@pytest.fixture
def sample_idx(tmp_path):
  # create a sample FASTA and index
  fasta = tmp_path / "sample.fasta"
  records = [
    SeqRecord(Seq("AAAAAA"), id="s1", description="s1"),
    SeqRecord(Seq("CCCCCC"), id="s2", description="s2"),
    SeqRecord(Seq("GGGGGG"), id="s3", description="s3"),
    SeqRecord(Seq("TTTTTTTT"), id="s4", description="s4"),  # Different length
  ]
  with open(fasta, "w") as fh:
    SeqIO.write(records, fh, "fasta")
  idx = tmp_path / "sample.idx"
  SeqIO.index_db(str(idx), str(fasta), "fasta")
  return str(idx)

def test_dataset_initialization(sample_idx):
  """Test valid Dataset initialization"""
  # Test DNA mode initialization
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)
  assert ds.mode == "dna"
  assert ds.kmer == 3
  assert ds.test_ratio == 0.5
  assert ds._tokenizer.vocab_size == 4
  assert ds.continuous == False  # default value

  # Test protein mode initialization
  ds_protein = Dataset(mode="protein", kmer=2, index_path=sample_idx)
  assert ds_protein.mode == "protein"
  assert ds_protein._tokenizer.vocab_size == 20
  assert ds_protein.test_ratio == 0.25  # default value

  # Test continuous mode
  ds_continuous = Dataset(mode="dna", kmer=3, index_path=sample_idx, continuous=True)
  assert ds_continuous.continuous == True
  assert ds_continuous._tokenizer.continuous == True

def test_initialization_errors(sample_idx):
  """Test Dataset initialization error handling"""
  # Missing kmer value
  with pytest.raises(ValueError, match="Must provide a kmer value"):
    Dataset(mode="dna", kmer=0, index_path=sample_idx)

  # Protein kmer size limit
  with pytest.raises(IndexError, match="Protein kmer size supported only up to 4"):
    Dataset(mode="protein", kmer=5, index_path=sample_idx)

  # DNA kmer size limit
  with pytest.raises(IndexError, match="DNA kmer size supported only up to 8"):
    Dataset(mode="dna", kmer=9, index_path=sample_idx)

  # Invalid file path
  with pytest.raises(FileNotFoundError):
    Dataset(mode="dna", kmer=3, index_path="nonexistent.idx")

def test_search_and_align(sample_idx):
  """Test sequence searching and alignment functionality"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)
  
  # Test search function
  groups = ds.search()
  assert isinstance(groups, dict)
  assert 6 in groups and 8 in groups
  assert set(groups[6]) == {"s1", "s2", "s3"}
  assert set(groups[8]) == {"s4"}

  # Test align function
  aligned = ds.align()
  assert isinstance(aligned, dict)
  assert aligned[6] == ["AAAAAA", "CCCCCC", "GGGGGG"]
  assert aligned[8] == ["TTTTTTTT"]

def test_fetch_sequence(sample_idx):
  """Test sequence fetching with sliding windows"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx)

  # Test overlapping windows
  windows = list(ds.fetch_sequence("s1", block_size=4, step=2))
  assert windows == ["AAAA", "AAAA", "AAAA"]  # overlapping windows from "AAAAAA"

  # Test non-overlapping windows (step defaults to block_size)
  windows = list(ds.fetch_sequence("s1", block_size=3))
  assert windows == ["AAA", "AAA"]

  # Test with different sequence
  windows = list(ds.fetch_sequence("s4", block_size=4, step=4))
  assert windows == ["TTTT", "TTTT"]  # non-overlapping from "TTTTTTTT"

  # Test error cases
  with pytest.raises(KeyError):
    list(ds.fetch_sequence("nonexistent", block_size=4))

  with pytest.raises(ValueError, match="shorter than block_size"):
    list(ds.fetch_sequence("s1", block_size=10))  # block_size > sequence length

def test_train_test_split(sample_idx):
  """Test train/test split functionality"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.25)
  train, val = ds.train_test_split()

  # With 4 items and ratio 0.25, train=3, val=1
  assert len(train) == 3 and len(val) == 1
  assert len(set(train + val)) == 4  # no overlap, all sequences included
  assert set(train + val) == {"s1", "s2", "s3", "s4"}

  # Test with different ratio
  ds_half = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)
  train_half, val_half = ds_half.train_test_split()
  assert len(train_half) == 2 and len(val_half) == 2

def test_get_batch_token_ids(sample_idx):
  """Test get_batch returning token IDs (one_hot=False)"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.0)
  batch = ds.get_batch(split="train", batch_size=2, block_size=4, one_hot=False)

  assert isinstance(batch, torch.Tensor)
  assert batch.shape == (2, 4)  # [batch_size, block_size]
  assert batch.dtype == torch.long

  # All values should be valid token IDs (0-3 for DNA)
  assert torch.all(batch >= 0) and torch.all(batch < 4)

def test_get_batch_one_hot(sample_idx):
  """Test get_batch returning one-hot encoded tensors (one_hot=True)"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.0)
  batch = ds.get_batch(split="train", batch_size=2, block_size=4, one_hot=True)

  assert isinstance(batch, torch.Tensor)
  assert batch.shape == (2, 4, 4)  # [batch_size, block_size, vocab_size]

  # One-hot: sum over vocab dimension should equal 1
  assert torch.all(batch.sum(-1) == 1)

  # All values should be 0 or 1
  assert torch.all((batch == 0) | (batch == 1))

def test_get_batch_validation_split(sample_idx):
  """Test get_batch with validation split"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)

  train_batch = ds.get_batch(split="train", batch_size=1, block_size=3, one_hot=False)
  val_batch = ds.get_batch(split="val", batch_size=1, block_size=3, one_hot=False)

  assert train_batch.shape == (1, 3)
  assert val_batch.shape == (1, 3)

def test_get_batch_errors(sample_idx):
  """Test get_batch error handling"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=1.0)  # All data in test

  # Should raise error when no training data available
  with pytest.raises(ValueError, match="No sequences available for split"):
    ds.get_batch(split="train", batch_size=1, block_size=3)

def test_get_sequence_stats(sample_idx):
  """Test sequence statistics"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx)
  stats = ds.get_sequence_stats()
  
  expected_keys = {"total_sequences", "min_length", "max_length", "mean_length", 
                   "vocab_size", "mode", "kmer", "continuous"}
  assert set(stats.keys()) == expected_keys

  assert stats["total_sequences"] == 4
  assert stats["min_length"] == 6
  assert stats["max_length"] == 8
  assert stats["mean_length"] == 6.5  # (6+6+6+8)/4
  assert stats["mode"] == "dna"
  assert stats["kmer"] == 3
  assert stats["vocab_size"] == 4
  assert stats["continuous"] == False

def test_dataset_length_and_str(sample_idx):
  """Test __len__ and __str__ methods"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx)

  assert len(ds) == 4
  str_repr = str(ds)
  assert "biosaic.Dataset" in str_repr
  assert "mode=dna" in str_repr
  assert "kmer=3" in str_repr
  assert "sequences=4" in str_repr

def test_continuous_mode(sample_idx):
  """Test continuous tokenization mode"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, continuous=True)
  assert ds.continuous == True
  assert ds._tokenizer.continuous == True

def test_protein_mode(sample_idx):
  """Test protein mode functionality"""
  ds = Dataset(mode="protein", kmer=2, index_path=sample_idx)
  assert ds.mode == "protein"
  assert ds._tokenizer.vocab_size == 20

  # Test that it can still process the sequences (even though they're DNA)
  # The tokenizer will just map unknown characters to 0
  batch = ds.get_batch(split="train", batch_size=1, block_size=3, one_hot=True)
  assert batch.shape == (1, 3, 20)  # protein vocab size is 20

def test_device_handling(sample_idx):
  """Test device handling in get_batch"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.0)
  
  # Test with CPU device (default)
  batch_cpu = ds.get_batch(split="train", batch_size=1, block_size=3, 
                          device=torch.device("cpu"))
  assert batch_cpu.device.type == "cpu"

def test_edge_cases(sample_idx):
  """Test edge cases and boundary conditions"""
  ds = Dataset(mode="dna", kmer=1, index_path=sample_idx)
  
  # Test with minimum kmer size
  assert ds.kmer == 1
  
  # Test with minimum block size
  batch = ds.get_batch(split="train", batch_size=1, block_size=1, one_hot=False)
  assert batch.shape == (1, 1)

# Self-executing test runner
if __name__ == "__main__":
  import tempfile
  import shutil
  from pathlib import Path
  
  print("Running Dataset tests...")
  
  # Create temporary directory for tests
  temp_dir = tempfile.mkdtemp()
  temp_path = Path(temp_dir)
  
  try:
    # Create sample index for manual testing
    fasta_path = temp_path / "test.fasta"
    records = [
      SeqRecord(Seq("AAAAAA"), id="seq1", description="Test sequence 1"),
      SeqRecord(Seq("CCCCCC"), id="seq2", description="Test sequence 2"),
      SeqRecord(Seq("GGGGGG"), id="seq3", description="Test sequence 3"),
      SeqRecord(Seq("TTTTTTTT"), id="seq4", description="Test sequence 4"),
    ]
    
    with open(fasta_path, "w") as fh:
      SeqIO.write(records, fh, "fasta")
    
    idx_path = temp_path / "test.idx"
    SeqIO.index_db(str(idx_path), str(fasta_path), "fasta")
    
    # Run basic functionality tests
    print("  ✓ Testing Dataset initialization...")
    ds = Dataset(mode="dna", kmer=3, index_path=str(idx_path))
    assert len(ds) == 4
    
    print("  ✓ Testing search and align...")
    groups = ds.search()
    aligned = ds.align()
    assert len(groups) == 2  # Two different lengths
    
    print("  ✓ Testing sequence fetching...")
    windows = list(ds.fetch_sequence("seq1", block_size=3))
    assert len(windows) == 2
    
    print("  ✓ Testing train/test split...")
    train, val = ds.train_test_split()
    assert len(train) + len(val) == 4
    
    print("  ✓ Testing batch generation...")
    batch = ds.get_batch(split="train", batch_size=2, block_size=3, one_hot=False)
    assert batch.shape == (2, 3)
    
    print("  ✓ Testing one-hot encoding...")
    batch_onehot = ds.get_batch(split="train", batch_size=2, block_size=3, one_hot=True)
    assert batch_onehot.shape == (2, 3, 4)
    
    print("  ✓ Testing statistics...")
    stats = ds.get_sequence_stats()
    assert stats["total_sequences"] == 4
    
    print("  ✓ Testing error handling...")
    try:
      Dataset(mode="dna", kmer=0, index_path=str(idx_path))
      assert False, "Should have raised ValueError"
    except ValueError:
      pass  # Expected

    try:
      Dataset(mode="protein", kmer=5, index_path=str(idx_path))
      assert False, "Should have raised IndexError"
    except IndexError:
      pass  # Expected

    try:
      Dataset(mode="dna", kmer=9, index_path=str(idx_path))
      assert False, "Should have raised IndexError"
    except IndexError:
      pass  # Expected

    try:
      # Use a path that doesn't exist but has the right format
      nonexistent_path = temp_path / "nonexistent.idx"
      Dataset(mode="dna", kmer=3, index_path=str(nonexistent_path))
      assert False, "Should have raised FileNotFoundError"
    except (FileNotFoundError, ValueError):
      pass  # Expected - could be either depending on BioPython version

    print("✓ All Dataset tests passed!")
  
  except Exception as e:
    print(f"  ✗ Test failed: {e}")
    import traceback
    traceback.print_exc()

  finally:
    # Close any open index files before cleanup
    try:
      if 'ds' in locals():
        ds._index.close()
      if 'ds_protein' in locals():
        ds_protein._index.close()
      if 'ds_continuous' in locals():
        ds_continuous._index.close()
    except:
      pass  # Ignore errors during cleanup
    
    # Cleanup with retry for Windows
    try:
      shutil.rmtree(temp_dir)
    except PermissionError:
      # Wait a bit and try again (Windows file locking issue)
      import time
      time.sleep(0.1)
      try:
        shutil.rmtree(temp_dir)
      except:
        print("  ⚠ Warning: Could not clean up temporary directory")
        pass
    print("  ✓ Cleanup completed.")