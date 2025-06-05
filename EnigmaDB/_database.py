import os
import time
import io
import sqlite3
import pandas as pd
from typing import List, Optional
from Bio import Entrez, SeqIO
from http.client import IncompleteRead
from urllib.error import HTTPError
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Database:
  """
  Automate DNA dataset collection & processing.
  Fetch raw FASTA from NCBI with guaranteed API key, immediate write
  and on-disk SeqIO index for O(1) lookup.
  
  Args:
    topics (List[str]): List containing strings of topics needed for building database
    out_dir (str): path to output directory.
    email (str): email address required by NCBI.
    api_key (str): NCBI API key for higher rate limits.
    max_rate (float): Max requests per second. Defaults to 3.0 (safe without API key).
    batch_size (int): Number of sequences per file. Defaults to 500.
    retmax (int): Max number of records to retrieve. Defaults to 10000.
    db (str): NCBI database to search. Defaults to 'nucleotide'.
    fmt (str): fetching format for the Search. defaults to 'fasta'.
  """

  def __init__(self, topics: List[str], out_dir: str, email: str, api_key: Optional[str] = None, 
               max_rate: float = 3.0, batch_size: int = 500, retmax: int = 10000, 
               db: str = 'nucleotide', fmt: str = 'fasta'):
    
    # Validate required parameters
    if not topics:
      raise ValueError("Topics list cannot be empty.")
    if not email:
      raise ValueError("An email address is required by NCBI.")
    if not out_dir:
      raise ValueError("Output directory must be specified.")
    
    # Set Entrez parameters
    Entrez.email = email
    if api_key:
      Entrez.api_key = api_key
      logger.info("Using NCBI API key for enhanced rate limits")
    else:
      logger.warning("No API key provided. Using default rate limits (max 3 requests/second)")
      max_rate = min(max_rate, 3.0)  # Enforce NCBI rate limit without API key

    # Validate parameters
    if max_rate <= 0:
      raise ValueError("max_rate must be positive")
    if batch_size <= 0:
      raise ValueError("batch_size must be positive")
    if retmax <= 0:
      raise ValueError("retmax must be positive")

    self.topics = topics
    self.out_dir = out_dir
    self.max_rate = max_rate
    self.batch_size = batch_size
    self.retmax = retmax
    self.db = db
    self.fmt = fmt
    self._sleep = 1.0 / self.max_rate
    
    # Create output directory
    try:
      os.makedirs(self.out_dir, exist_ok=True)
      logger.info(f"Output directory: {self.out_dir}")
    except OSError as e:
      raise OSError(f"Failed to create output directory: {e}")

  def _sanitize(self, s: str) -> str:
    """Sanitize filename by removing invalid characters."""
    if not s:
      return "unnamed"
    sanitized = "".join(c if c.isalnum() or c in "-_" else "_" for c in s).strip("_")
    return sanitized if sanitized else "unnamed"

  def search(self, query: str) -> List[str]:
    """Search NCBI database for sequence IDs matching the query."""
    if not query.strip():
      logger.warning("Empty query provided")
      return []
    
    try:
      logger.info(f"Searching for: {query}")
      handle = Entrez.esearch(db=self.db, term=query, retmax=self.retmax)
      rec = Entrez.read(handle)
      handle.close()
      
      id_list = rec.get('IdList', [])
      logger.info(f"Found {len(id_list)} sequences for query: {query}")
      return id_list
      
    except Exception as e:
      logger.error(f"Search failed for query '{query}': {e}")
      return []

  def _safe_efetch(self, ids: List[str]) -> io.StringIO:
    """Fetch sequences with retry logic."""
    if not ids:
      raise ValueError("No IDs provided for fetching")
    
    for attempt in range(3):
      try:
        logger.debug(f"Fetching batch of {len(ids)} sequences (attempt {attempt + 1})")
        handle = Entrez.efetch(db=self.db, id=",".join(ids), rettype=self.fmt, retmode='text')
        data = handle.read()
        handle.close()
        
        if not data:
          raise ValueError("Empty response from NCBI")
        
        return io.StringIO(data)
        
      except (IncompleteRead, HTTPError, Exception) as e:
        logger.warning(f"Fetch attempt {attempt + 1} failed: {e}")
        if attempt < 2:  # Don't sleep on the last attempt
          time.sleep(2 ** attempt)
    
    raise RuntimeError(f"Failed to fetch batch of {len(ids)} sequences after 3 attempts")

  def build(self, with_index: bool = True):
    """
    For each topic:
      1. search -> get UIDs
      2. fetch in batches, write as you go
      3. build an on-disk index for O(1) lookups
    """
    successful_builds = 0
    
    for i, topic in enumerate(self.topics, 1):
      logger.info(f"[{i}/{len(self.topics)}] Processing topic: {topic}")
      
      try:
        # Search for sequences
        ids = self.search(topic)
        if not ids:
          logger.warning(f"No sequences found for topic: {topic}")
          continue

        # Prepare file paths
        fname = self._sanitize(topic)
        fasta_path = os.path.join(self.out_dir, f"{fname}.fasta")
        idx_path = os.path.join(self.out_dir, f"{fname}.idx")

        # Check if files already exist
        if os.path.exists(fasta_path):
          logger.info(f"FASTA file already exists: {fasta_path}")
          user_input = input(f"Overwrite {fasta_path}? (y/n): ").lower().strip()
          if user_input != 'y':
            logger.info("Skipping existing file")
            continue

        sequences_written = 0
        
        # Fetch sequences in batches
        with open(fasta_path, 'w', encoding='utf-8') as out_fh:
          for batch_start in range(0, len(ids), self.batch_size):
            batch_end = min(batch_start + self.batch_size, len(ids))
            batch = ids[batch_start:batch_end]
            
            try:
              logger.info(f"Fetching batch {batch_start//self.batch_size + 1}/{(len(ids)-1)//self.batch_size + 1}")
              handle = self._safe_efetch(batch)
              
              # Count sequences in this batch
              batch_content = handle.getvalue()
              batch_seq_count = batch_content.count('>')
              
              # Write batch to file
              handle.seek(0)  # Reset StringIO position
              out_fh.write(batch_content)
              out_fh.flush()
              
              sequences_written += batch_seq_count
              logger.debug(f"Wrote {batch_seq_count} sequences from batch")
              
            except Exception as e:
              logger.error(f"Batch {batch_start}-{batch_end} failed: {e}")
              continue
            
            # Rate limiting
            time.sleep(self._sleep)

        logger.info(f"Wrote {sequences_written} sequences to {fasta_path}")

        # Build index if requested
        if with_index and os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
          try:
            logger.info(f"Building index: {idx_path}")
            SeqIO.index_db(idx_path, fasta_path, "fasta")
            logger.info(f"Successfully built index: {idx_path}")
          except Exception as e:
            logger.error(f"Could not build index for {fasta_path}: {e}")

        successful_builds += 1
        logger.info(f"Completed processing topic: {topic}\n")
        
      except Exception as e:
        logger.error(f"Failed to process topic '{topic}': {e}")
        continue

    logger.info(f"Successfully processed {successful_builds}/{len(self.topics)} topics")

def create_index(input_dir: str, index_path: str = "combined.idx"):
  """
  Build a single O(1) lookup index (.idx) from all FASTA files in a directory,
  skipping duplicate record IDs automatically.

  Args:
    input_dir: Directory containing .fasta/.fa files.
    index_path: Path to SQLite index file to create.
  """
  if not os.path.exists(input_dir):
    raise FileNotFoundError(f"Input directory not found: {input_dir}")
  
  if not os.path.isdir(input_dir):
    raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

  # Gather FASTA files
  fasta_files = []
  for fn in os.listdir(input_dir):
    if fn.lower().endswith(('.fa', '.fasta')):
      full_path = os.path.join(input_dir, fn)
      if os.path.isfile(full_path) and os.path.getsize(full_path) > 0:
        fasta_files.append(full_path)

  if not fasta_files:
    logger.warning(f"No FASTA files found in '{input_dir}'")
    return

  logger.info(f"Found {len(fasta_files)} FASTA files to index")

  # Create SQLite DB and table
  try:
    # Remove existing index if it exists
    if os.path.exists(index_path):
      os.remove(index_path)
    
    conn = sqlite3.connect(index_path)
    c = conn.cursor()
    c.execute("""
      CREATE TABLE seq_index (
        key      TEXT PRIMARY KEY,
        filename TEXT NOT NULL,
        offset   INTEGER NOT NULL,
        length   INTEGER NOT NULL
      )
    """)
    conn.commit()

    total_sequences = 0
    duplicate_count = 0

    # Scan each FASTA and record offsets
    for fasta_file in fasta_files:
      logger.info(f"Indexing: {os.path.basename(fasta_file)}")
      
      try:
        with open(fasta_file, 'r', encoding='utf-8') as fh:
          while True:
            pos = fh.tell()
            header = fh.readline()
            
            if not header:
              break
            
            if not header.startswith('>'):
              continue
            
            # Extract sequence ID (first word after >)
            seq_id = header[1:].split()[0].strip()
            if not seq_id:
              logger.warning(f"Empty sequence ID at position {pos} in {fasta_file}")
              continue
            
            # Read sequence lines until next '>' or EOF
            seq_start_pos = fh.tell()
            seq_len = 0
            
            while True:
              line_start = fh.tell()
              line = fh.readline()
              
              if not line or line.startswith('>'):
                # Rewind if we overshot onto next header
                if line and line.startswith('>'):
                  fh.seek(line_start)
                break
              
              seq_len += len(line.strip())
            
            # Insert or ignore duplicates
            try:
              c.execute("INSERT INTO seq_index (key, filename, offset, length) VALUES (?, ?, ?, ?)", 
                       (seq_id, fasta_file, pos, seq_len))
              total_sequences += 1
            except sqlite3.IntegrityError:
              duplicate_count += 1
              logger.debug(f"Duplicate sequence ID ignored: {seq_id}")

      except Exception as e:
        logger.error(f"Error processing {fasta_file}: {e}")
        continue

    conn.commit()
    conn.close()
    
    logger.info(f"Index creation complete:")
    logger.info(f"  - Total sequences indexed: {total_sequences}")
    logger.info(f"  - Duplicates skipped: {duplicate_count}")
    logger.info(f"  - Index saved to: {index_path}")
    
  except Exception as e:
    logger.error(f"Failed to create index: {e}")
    if 'conn' in locals():
      conn.close()
    if os.path.exists(index_path):
      os.remove(index_path)
    raise

def convert_fasta(input_dir: str, output_dir: str, mode: str = 'csv'):
  """
  Read all FASTA files in `input_dir` and write out CSV or Parquet files
  with columns: id, name (description), length (original), sequence.

  Args:
    input_dir: Path to folder containing .fasta files.
    output_dir: Path where output files will be saved.
    mode: 'csv' or 'parquet'.
  """
  if not os.path.exists(input_dir):
    raise FileNotFoundError(f"Input directory not found: {input_dir}")
  
  if mode not in ['csv', 'parquet']:
    raise ValueError(f"Unsupported mode: {mode}. Use 'csv' or 'parquet'")

  try:
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
  except OSError as e:
    raise OSError(f"Failed to create output directory: {e}")

  # Check for required dependencies
  if mode == 'parquet':
    try:
      import pyarrow
    except ImportError:
      raise ImportError("pyarrow is required for parquet output. Install with: pip install pyarrow")

  converted_files = 0
  
  for fname in os.listdir(input_dir):
    if not fname.lower().endswith(('.fasta', '.fa')):
      continue
    
    path_in = os.path.join(input_dir, fname)
    
    if not os.path.isfile(path_in) or os.path.getsize(path_in) == 0:
      logger.warning(f"Skipping empty or invalid file: {fname}")
      continue
    
    try:
      logger.info(f"Converting: {fname}")
      
      # Parse FASTA file
      records = list(SeqIO.parse(path_in, 'fasta'))
      
      if not records:
        logger.warning(f"No sequences found in: {fname}")
        continue
      
      # converting to DataFrame
      rows = []
      for record in records:
        seq_str = str(record.seq)
        rows.append({
          'id': record.id,
          'name': record.description,
          'length': len(seq_str),
          'sequence': seq_str
        })

      df = pd.DataFrame(rows)
      base = os.path.splitext(fname)[0] # generating output path
      
      if mode == 'csv':
        out_path = os.path.join(output_dir, f"{base}.csv")
        df.to_csv(out_path, index=False)
      elif mode == 'parquet':
        out_path = os.path.join(output_dir, f"{base}.parquet")
        df.to_parquet(out_path, index=False)
      
      logger.info(f"Converted {len(records)} sequences: {path_in} -> {out_path}")
      converted_files += 1
      
    except Exception as e:
      logger.error(f"Error converting {fname}: {e}")
      continue

  logger.info(f"Conversion complete. Processed {converted_files} files")