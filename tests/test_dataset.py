import torch
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from EnigmaDB import Dataset

# stub tokenizer
class DummyTokenizer:
  vocab_size = 4
  def encode(self, s):
    # map A,C,G,T → 0,1,2,3
    m = {'A':0,'C':1,'G':2,'T':3,'N':0}
    return [m.get(ch,0) for ch in s]

@pytest.fixture(autouse=True)
def patch_tokenizer(monkeypatch):
  # monkey-patch biosaic.tokenizer and get_encodings
  import builtins
  monkeypatch.setitem(__import__('sys').modules, 'biosaic', type('m',(),{})())
  import biosaic
  biosaic.tokenizer = lambda encoding: DummyTokenizer()
  biosaic.get_encodings = [None]*6
  yield

@pytest.fixture
def sample_idx(tmp_path):
  # create a sample FASTA and index
  fasta = tmp_path / "sample.fasta"
  records = [
    SeqRecord(Seq("AAAAAA"), id="s1", description="s1"),
    SeqRecord(Seq("CCCCCC"), id="s2", description="s2"),
    SeqRecord(Seq("GGGGGG"), id="s3", description="s3"),
  ]
  with open(fasta, "w") as fh:
    SeqIO.write(records, fh, "fasta")
  idx = tmp_path / "sample.idx"
  SeqIO.index_db(str(idx), str(fasta), "fasta")
  return str(idx)

def test_search_and_align(sample_idx):
  ds = Dataset(kmer=3, index_path=sample_idx, test_ratio=0.5)
  groups = ds.search()
  # all sequences length 6
  assert list(groups.keys()) == [6]
  assert set(groups[6]) == {"s1","s2","s3"}
  aligned = ds.align()
  assert aligned[6] == ["AAAAAA","CCCCCC","GGGGGG"]

def test_fetch_sequence(sample_idx):
  ds = Dataset(kmer=3, index_path=sample_idx)
  windows = list(ds.fetch_sequence("s1", block_size=4, step=2))
  assert windows == ["AAAA","AAAA","AAAA"]  # overlapping windows

def test_train_test_split(sample_idx):
  ds = Dataset(kmer=3, index_path=sample_idx, test_ratio=0.33)
  train, val = ds.train_test_split()
  # with 3 items and ratio 0.33, train=2, val=1
  assert len(train) == 2 and len(val) == 1

def test_get_batch_shapes(sample_idx):
  ds = Dataset(kmer=3, index_path=sample_idx, test_ratio=0.0)
  batch = ds.get_batch(split="train", batch_size=2, block_size=4)
  assert isinstance(batch, torch.Tensor)
  assert batch.shape == (2,4,ds._tokenizer.vocab_size)
  # values one-hot: sum over vocab dim ==1
  assert torch.all(batch.sum(-1)==1)