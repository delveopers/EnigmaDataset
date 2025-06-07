import torch
import pytest, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from EnigmaDB import Dataset

class DummyTokenizer:
  def __init__(self, mode="dna", kmer=1, continuous=False):
    self.mode = mode
    self.kmer = kmer
    self.continuous = continuous
    self.vocab_size = 4 if mode == "dna" else 20  # DNA: 4 bases, Protein: 20 amino acids
    
  def encode(self, s):
    # map A,C,G,T → 0,1,2,3 for DNA
    if self.mode == "dna":
      m = {'A':0,'C':1,'G':2,'T':3,'N':0}
    else:  # protein mode
      # Simple mapping for amino acids
      amino_acids = "ACDEFGHIKLMNPQRSTVWY"
      m = {aa: i for i, aa in enumerate(amino_acids)}
    return [m.get(ch.upper(), 0) for ch in s]

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
  # Test valid initialization
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)
  assert ds.mode == "dna"
  assert ds.kmer == 3
  assert ds.test_ratio == 0.5
  assert ds._tokenizer.vocab_size == 4

  # Test protein mode initialization
  ds_protein = Dataset(mode="protein", kmer=2, index_path=sample_idx)
  assert ds_protein.mode == "protein"
  assert ds_protein._tokenizer.vocab_size == 20

def test_initialization_errors(sample_idx):
  # Test kmer validation
  with pytest.raises(ValueError, match="Must provide a kmer value"):
    Dataset(mode="dna", kmer=0, index_path=sample_idx)

  # Test protein kmer size limit
  with pytest.raises(IndexError, match="Protein kmer size supported only up to 4"):
    Dataset(mode="protein", kmer=5, index_path=sample_idx)

  # Test DNA kmer size limit
  with pytest.raises(IndexError, match="DNA kmer size supported only up to 8"):
    Dataset(mode="dna", kmer=9, index_path=sample_idx)

  # Test invalid file path
  with pytest.raises(FileNotFoundError):
    Dataset(mode="dna", kmer=3, index_path="nonexistent.idx")

def test_search_and_align(sample_idx):
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.5)
  groups = ds.search()
  # sequences of length 6 and 8
  assert 6 in groups and 8 in groups
  assert set(groups[6]) == {"s1", "s2", "s3"}
  assert set(groups[8]) == {"s4"}

  aligned = ds.align()
  assert aligned[6] == ["AAAAAA", "CCCCCC", "GGGGGG"]
  assert aligned[8] == ["TTTTTTTT"]

def test_fetch_sequence(sample_idx):
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx)

  # Test overlapping windows
  windows = list(ds.fetch_sequence("s1", block_size=4, step=2))
  assert windows == ["AAAA", "AAAA", "AAAA"]  # overlapping windows

  # Test non-overlapping windows
  windows = list(ds.fetch_sequence("s1", block_size=3))  # step defaults to block_size
  assert windows == ["AAA", "AAA"]

  # Test error cases
  with pytest.raises(KeyError):
      list(ds.fetch_sequence("nonexistent", block_size=4))

  with pytest.raises(ValueError):
    list(ds.fetch_sequence("s1", block_size=10))  # block_size > sequence length

def test_train_test_split(sample_idx):
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx, test_ratio=0.25)
  train, val = ds.train_test_split()

  # with 4 items and ratio 0.25, train=3, val=1
  assert len(train) == 3 and len(val) == 1
  assert len(set(train + val)) == 4  # no overlap

  # Test empty index error would be raised by SeqIO, not our code
  # but we can test the behavior with our current setup

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
  with pytest.raises(ValueError):
    ds.get_batch(split="train", batch_size=1, block_size=3)

def test_get_sequence_stats(sample_idx):
  """Test sequence statistics"""
  ds = Dataset(mode="dna", kmer=3, index_path=sample_idx)
  stats = ds.get_sequence_stats()
  expected_keys = {"total_sequences", "min_length", "max_length", "mean_length", "vocab_size", "mode", "kmer", "continuous"}
  assert set(stats.keys()) == expected_keys

  assert stats["total_sequences"] == 4
  assert stats["min_length"] == 6
  assert stats["max_length"] == 8
  assert stats["mode"] == "dna"
  assert stats["kmer"] == 3
  assert stats["vocab_size"] == 4

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