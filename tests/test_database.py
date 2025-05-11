import os
import io
import sqlite3
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from EnigmaDB import Database, create_index, convert_fasta

class DummyDatabase(Database):
  """ Override network methods for testing. """
  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def search(self, query):
    # pretend we found three fake IDs
    return ['id1','id2','id3']

  def _safe_efetch(self, ids):
    # return a small FASTA stream
    records = [
      SeqRecord(Seq("AAA"), id=ids[0], description="desc1"),
      SeqRecord(Seq("CCC"), id=ids[1], description="desc2"),
    ]
    stream = io.StringIO()
    SeqIO.write(records, stream, "fasta")
    stream.seek(0)
    return stream

@pytest.fixture
def tmp_out(tmp_path):
  d = tmp_path / "out"
  d.mkdir()
  return str(d)

def test_sanitize(tmp_out):
  db = DummyDatabase(['Q'], tmp_out, email="a@b.com", api_key="KEY")
  assert db._sanitize("Hello World!") == "Hello_World_"
  assert db._sanitize("A/B:C") == "A_B_C"

def test_build_and_index(tmp_out):
  db = DummyDatabase(['TestTopic'], tmp_out, email="a@b.com", api_key="KEY", batch_size=2)
  db.build(with_index=True)
  fasta = os.path.join(tmp_out, "TestTopic.fasta")
  idx   = os.path.join(tmp_out, "TestTopic.idx")
  # FASTA file exists and contains our dummy records
  assert os.path.exists(fasta)
  recs = list(SeqIO.parse(fasta, "fasta"))
  assert {r.id for r in recs} == {'id1','id2'}
  # idx file should be a SQLite DB
  conn = sqlite3.connect(idx)
  cur = conn.cursor()
  cur.execute("SELECT key,filename,offset,length FROM seq_index")
  rows = cur.fetchall()
  assert len(rows) == 2
  keys = {r[0] for r in rows}
  assert keys == {'id1','id2'}
  conn.close()

def test_create_index_and_convert(tmp_path):
  # create two tiny FASTAs
  dir_in = tmp_path / "fa"
  dir_in.mkdir()
  for name,seq in [("a", "ATG"), ("b","GGA")]:
    path = dir_in / f"{name}.fasta"
    rec = SeqRecord(Seq(seq), id=name, description=name)
    with open(path, "w") as fh:
      SeqIO.write([rec], fh, "fasta")
  idx_file = tmp_path / "combined.idx"
  # build combined index
  create_index(str(dir_in), str(idx_file))
  assert idx_file.exists()
  # test lookup from index
  from Bio import SeqIO as sio
  index = sio.index_db(str(idx_file))
  assert "a" in index and "b" in index
  assert str(index["a"].seq) == "ATG"
  index.close()
  # convert to CSV
  out_dir = tmp_path / "out"
  convert_fasta(str(dir_in), str(out_dir), mode="csv")
  csv_files = list(out_dir.glob("*.csv"))
  assert len(csv_files) == 2
  df = pd.read_csv(csv_files[0])
  assert set(df.columns) == {"id","name","length","sequence"}