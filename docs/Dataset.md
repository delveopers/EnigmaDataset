# 🧬 `Dataset` Class — Biosaic Batch Loader

The `Dataset` class is a specialized data loader for training machine learning models on biological sequences (DNA/Protein). It provides efficient batch generation with k-mer tokenization and supports both DNA and protein sequences through the biosaic tokenizer system.

## Features

* **Efficient O(1) indexed lookups** via Biopython's `SeqIO.index_db`
* **K-mer tokenization** for DNA (up to 8-mer) and protein (up to 4-mer) sequences
* **Flexible batch generation** with token IDs or one-hot encoding
* **Sliding window sequence extraction** with configurable step sizes
* **Automatic train/test splitting** with configurable ratios
* **Continuous tokenization support** for overlapping k-mers

## Class Signature

```python
Dataset(mode: str, kmer: int, index_path: str, continuous: bool = False, test_ratio: float = 0.25)
```

### Parameters

| Parameter    | Type    | Description                                           |
|------------- |---------|-------------------------------------------------------|
| `mode`       | `str`   | Sequence type: `"dna"` or `"protein"`                |
| `kmer`       | `int`   | K-mer size (DNA: ≤8, Protein: ≤4)                   |
| `index_path` | `str`   | Path to FASTA file or existing `.idx` file           |
| `continuous` | `bool`  | Enable sliding window tokenization (default: False)  |
| `test_ratio` | `float` | Proportion for test split (default: 0.25)            |

## Dependencies

* [Biopython](https://biopython.org/) - For FASTA file handling
* [NumPy](https://numpy.org/) - For numerical operations
* [PyTorch](https://pytorch.org/) - For tensor operations
* [Biosaic](https://github.com/your-repo/biosaic) - Custom tokenizer system

## Core Methods

### `search() → Dict[int, List[str]]`

Groups sequence IDs by their exact length for efficient batch processing.

```python
>>> dataset.search()
{150: ['seq1', 'seq3'], 200: ['seq2', 'seq4']}
```

### `align() → Dict[int, List[str]]`

Returns raw sequence strings grouped by identical lengths.

```python
>>> dataset.align()
{150: ['ATGCGT...', 'CGTAAT...'], 200: ['ATGCGTAA...', 'CGTAATGG...']}
```

### `fetch_sequence(seq_id: str, block_size: int, step: int = None) → Iterator[str]`

Yields sliding windows from a sequence with specified block size and step.

```python
>>> list(dataset.fetch_sequence("seq1", block_size=10, step=5))
['ATGCGTAACG', 'GTAACGTCGA', 'ACGTCGATAG', ...]
```

**Parameters:**
- `seq_id`: Sequence identifier from the index
- `block_size`: Size of each window
- `step`: Step size between windows (default: `block_size` for non-overlapping)

### `train_test_split() → Tuple[List[str], List[str]]`

Randomly splits sequence IDs into training and testing sets.

```python
>>> train_ids, test_ids = dataset.train_test_split()
>>> len(train_ids), len(test_ids)
(750, 250)  # For 1000 sequences with test_ratio=0.25
```

### `get_batch(split: str, batch_size: int, block_size: int, step: int = None, device: torch.device = torch.device("cpu"), one_hot: bool = False) → torch.Tensor`

Generates batches of tokenized sequences for training or evaluation.

```python
# Token ID batch
batch = dataset.get_batch("train", batch_size=32, block_size=128)
# Shape: [32, 128] - token IDs

# One-hot encoded batch  
batch_onehot = dataset.get_batch("train", batch_size=32, block_size=128, one_hot=True)
# Shape: [32, 128, vocab_size] - one-hot vectors
```

**Parameters:**
- `split`: `"train"` or `"val"`/`"test"`
- `batch_size`: Number of sequences in the batch
- `block_size`: Length of each sequence window
- `step`: Step size for sliding windows (optional)
- `device`: Target device for tensors
- `one_hot`: Return one-hot encoded tensors if True

### `get_sequence_stats() → Dict[str, Any]`

Returns comprehensive statistics about the loaded dataset.

```python
>>> dataset.get_sequence_stats()
{
    'total_sequences': 1000,
    'min_length': 50,
    'max_length': 2000,
    'mean_length': 425.6,
    'vocab_size': 64,  # For 3-mer DNA
    'mode': 'dna',
    'kmer': 3,
    'continuous': False
}
```

## Example Usage

### Basic DNA Dataset

```python
from biosaic.dataset import Dataset

# Initialize dataset for DNA sequences
dataset = Dataset(
    mode="dna",
    kmer=3,
    index_path="data/sequences.fasta",
    continuous=True,
    test_ratio=0.2
)

# Get dataset statistics
stats = dataset.get_sequence_stats()
print(f"Dataset contains {stats['total_sequences']} sequences")
print(f"Vocabulary size: {stats['vocab_size']}")

# Generate training batch
train_batch = dataset.get_batch(
    split="train",
    batch_size=16,
    block_size=256,
    device=torch.device("cuda" if torch.cuda.is_available() else "cpu")
)
print(f"Batch shape: {train_batch.shape}")  # [16, 256]
```

### Protein Dataset with One-Hot Encoding

```python
# Initialize protein dataset
protein_dataset = Dataset(
    mode="protein",
    kmer=2,
    index_path="data/proteins.fasta"
)

# Get one-hot encoded batch for validation
val_batch = protein_dataset.get_batch(
    split="val",
    batch_size=8,
    block_size=100,
    one_hot=True
)
print(f"One-hot batch shape: {val_batch.shape}")  # [8, 100, vocab_size]
```

### Working with Sequence Windows

```python
# Explore sequence windows
for window in dataset.fetch_sequence("sequence_001", block_size=50, step=10):
    print(f"Window: {window[:20]}...")  # First 20 characters
    if len(list(dataset.fetch_sequence("sequence_001", 50, 10))) > 5:
        break  # Show first 5 windows only
```

## Error Handling

The class includes comprehensive error handling for common issues:

- **FileNotFoundError**: Invalid FASTA file path
- **ValueError**: Block size larger than sequence length, empty datasets
- **IndexError**: K-mer size exceeds supported limits
- **KeyError**: Invalid sequence ID requests

## Performance Considerations

* **Memory Efficiency**: Uses on-disk indexing for large datasets
* **I/O Optimization**: SeqIO.index_db provides O(1) sequence access
* **Batch Generation**: Efficient random sampling without loading entire dataset
* **Tokenizer Integration**: Direct integration with biosaic tokenizer for consistent encoding

## Notes

* DNA sequences support k-mer sizes up to 8 (vocab size: 4^k)
* Protein sequences support k-mer sizes up to 4 (vocab size: 20^k)
* The `continuous` parameter enables overlapping k-mer extraction
* All sequences are automatically converted to uppercase
* Index files are automatically created if they don't exist

## Building FASTA Indexes

To create an index file for your FASTA data:

```python
from Bio import SeqIO

# Create index from single FASTA file
SeqIO.index_db("sequences.idx", "sequences.fasta", "fasta")

# Create index from multiple FASTA files
SeqIO.index_db("combined.idx", ["file1.fasta", "file2.fasta"], "fasta")
```
