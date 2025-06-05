"""
  @dataset.py
    * contains Dataset class: special class for loading, formatting & creating batches of datasets
     * compatible with biosaic tokenizer system
    * Raises:
        FileNotFoundError: invalid file path is given
        ValueError: data length is less than block size
        ValueError: if data is not loaded for performing train*test split
        IndexError: out of range index
    * Returns:
        torch.tensor(): return batches of tokenized DNA/Protein datasets
"""

import torch
import numpy as np
from Bio import SeqIO
from typing import List, Dict, Iterator
from biosaic import Tokenizer

class Dataset:
  """
  Dataset backed by an on-disk FASTA index for O(1) sequence lookup.

  Functions:
    - search(): group sequence IDs by exact length -> Dict[int, List[str]]
    - align(): pad each group's sequences to its max length -> Dict[int, List[str]]
    - fetch_sequence(): sliding windows from a single sequence
    - train_test_split(): split IDs into train/test lists
    - get_batch(): produce encoded batches [B, block_size] or one-hot [B, block_size, vocab_size]
  """

  def __init__(self, mode: str, kmer: int, index_path: str, continuous: bool = False, test_ratio: float = 0.25):
    """
    Initialize Dataset with tokenizer and FASTA index.
    
    Args:
      mode (str): "dna" or "protein" - type of sequences
      kmer (int): k-mer size for tokenization
      index_path (str): path to FASTA file or index
      continuous (bool): whether to use sliding window tokenization
      test_ratio (float): fraction of data to use for testing
    """
    if not kmer:
      raise ValueError("Must provide a kmer value!")
    
    # Validate kmer limits based on mode
    if mode == "protein" and kmer > 4:
      raise IndexError("Protein kmer size supported only up to 4!")
    elif mode == "dna" and kmer > 8:
      raise IndexError("DNA kmer size supported only up to 8!")
    
    self.mode = mode
    self.kmer = kmer
    self.continuous = continuous
    self.test_ratio = test_ratio
    
    # Initialize tokenizer using the main Tokenizer class
    self._tokenizer = Tokenizer(mode=mode, kmer=kmer, continuous=continuous)
    
    # Build or load on-disk index
    try:
      self._index = SeqIO.index_db(index_path)
    except FileNotFoundError:
      raise FileNotFoundError(f"FASTA file not found at: {index_path}")

  def search(self) -> Dict[int, List[str]]:
    """
    Group sequence IDs by their exact length.
    Returns: { length: [seq_id, ...], ... }
    """
    groups: Dict[int, List[str]] = {}
    for seq_id in self._index.keys():
      length = len(self._index[seq_id].seq)
      groups.setdefault(length, []).append(seq_id)
    return groups

  def align(self) -> Dict[int, List[str]]:
    """
    For each length group, return raw sequence strings.
    Returns: { length: [sequence_str, ...], ... }
    """
    aligned: Dict[int, List[str]] = {}
    for length, ids in self.search().items():
      seqs = []
      for seq_id in ids:
        seqs.append(str(self._index[seq_id].seq).upper())
      aligned[length] = seqs
    return aligned

  def fetch_sequence(self, seq_id: str, block_size: int, step: int = None) -> Iterator[str]:
    """
    Yield sliding windows of size `block_size` from sequence `seq_id`.
    
    Args:
      seq_id: sequence identifier
      block_size: size of each window
      step: increment size; defaults to block_size (non-overlapping)
    """
    if seq_id not in self._index:
      raise KeyError(f"Sequence ID '{seq_id}' not found in index")
    
    seq = str(self._index[seq_id].seq).upper()
    N = len(seq)
    
    if block_size > N:
      raise ValueError(f"Sequence {seq_id} (length {N}) is shorter than block_size ({block_size})")
    
    if step is None:
      step = block_size
    
    for i in range(0, N - block_size + 1, step):
      yield seq[i:i + block_size]

  def train_test_split(self) -> tuple[List[str], List[str]]:
    """
    Split the full list of seq_ids into train/test according to test_ratio.
    Returns: (train_ids, test_ids)
    """
    all_ids = list(self._index.keys())
    if not all_ids:
      raise ValueError("No sequences found in the index")
    
    # Shuffle for random split
    np.random.shuffle(all_ids)
    split_idx = int(len(all_ids) * (1 - self.test_ratio))
    return all_ids[:split_idx], all_ids[split_idx:]

  def get_batch(self, split: str, batch_size: int, block_size: int, step: int = None, 
                device: torch.device = torch.device("cpu"), one_hot: bool = False) -> torch.Tensor:
    """
    Returns encoded batch tensor.
    
    Args:
      split: "train" or "val"/"test"
      batch_size: number of sequences per batch
      block_size: length of each sequence window
      step: step size for sliding windows
      device: torch device to place tensors on
      one_hot: if True, return one-hot encoded tensors [B, block_size, vocab_size]
               if False, return token IDs [B, block_size]
    
    Returns:
      torch.Tensor: batch of encoded sequences
    """
    # Get train/test split
    train_ids, val_ids = self.train_test_split()
    ids = train_ids if split == "train" else val_ids
    
    if not ids:
      raise ValueError(f"No sequences available for split: {split}")

    # Collect samples
    samples = []
    for _ in range(batch_size):
      # Randomly select a sequence
      seq_id = np.random.choice(ids)
      
      # Get all possible windows
      try:
        windows = list(self.fetch_sequence(seq_id, block_size, step))
        if not windows:
          raise ValueError(f"No valid windows for sequence {seq_id}")
        
        # Pick a random window
        subseq = np.random.choice(windows)
        
        # Tokenize the subsequence
        token_ids = self._tokenizer.encode(subseq)
        
        # Ensure we have the right length
        if len(token_ids) == 0:
          raise ValueError(f"Empty tokenization for sequence: {subseq}")
        
        # Convert to tensor
        token_tensor = torch.tensor(token_ids, dtype=torch.long)
        
        if one_hot:
          # Convert to one-hot encoding
          oh = torch.nn.functional.one_hot(
            token_tensor,
            num_classes=self._tokenizer.vocab_size
          )
          samples.append(oh)
        else:
          samples.append(token_tensor)
          
      except (ValueError, KeyError) as e:
        # Skip problematic sequences and try another
        continue

    if not samples:
      raise ValueError("Could not generate any valid samples")

    # Stack into batch tensor
    try:
      batch = torch.stack(samples).to(device)
      return batch
    except RuntimeError as e:
      raise ValueError(f"Error stacking batch tensors: {e}")

  def get_sequence_stats(self) -> Dict[str, any]:
    """
    Get statistics about the loaded sequences.
    
    Returns:
      Dict containing sequence statistics
    """
    if not hasattr(self, '_index') or not self._index:
      raise ValueError("No sequences loaded")
    
    lengths = [len(self._index[seq_id].seq) for seq_id in self._index.keys()]
    
    return {
      "total_sequences": len(lengths),
      "min_length": min(lengths) if lengths else 0,
      "max_length": max(lengths) if lengths else 0,
      "mean_length": np.mean(lengths) if lengths else 0,
      "vocab_size": self._tokenizer.vocab_size,
      "mode": self.mode,
      "kmer": self.kmer,
      "continuous": self.continuous
    }

  def __len__(self):
    """Return number of sequences in the dataset."""
    return len(self._index)

  def __str__(self):
    return f"biosaic.Dataset <mode={self.mode}, kmer={self.kmer}, continuous={self.continuous}, sequences={len(self)}>"