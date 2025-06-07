# 🧬 NCBI Database Builder

A comprehensive Python utility for **downloading**, **organizing**, and **processing** biological sequence datasets from NCBI. The `Database` class automates the entire pipeline from query construction to indexed FASTA file generation.

## Features

* **Automated NCBI Downloads**: Fetch sequences using Entrez API with rate limiting
* **Batch Processing**: Configurable batch sizes with retry logic for robust downloads
* **On-Disk Indexing**: Automatic SeqIO index creation for O(1) sequence access
* **Format Conversion**: Convert FASTA files to CSV or Parquet formats
* **Index Management**: Create combined indexes from multiple FASTA files
* **Comprehensive Logging**: Detailed progress tracking and error reporting

## Installation

```bash
pip install biopython pandas pyarrow sqlite3
```

## Core Classes and Functions

### `Database` Class

Main class for automated NCBI sequence retrieval and processing.

```python
Database(
    topics: List[str], 
    out_dir: str, 
    email: str, 
    api_key: Optional[str] = None,
    max_rate: float = 3.0, 
    batch_size: int = 500, 
    retmax: int = 10000,
    db: str = 'nucleotide', 
    fmt: str = 'fasta'
)
```

#### Parameters

| Parameter    | Type           | Description                                    | Default      |
|-------------|----------------|------------------------------------------------|--------------|
| `topics`    | `List[str]`    | List of Entrez search queries                  | *required*   |
| `out_dir`   | `str`          | Output directory for downloaded files          | *required*   |
| `email`     | `str`          | Email address (required by NCBI)              | *required*   |
| `api_key`   | `str`          | NCBI API key for enhanced rate limits          | `None`       |
| `max_rate`  | `float`        | Maximum requests per second                    | `3.0`        |
| `batch_size`| `int`          | Sequences per download batch                   | `500`        |
| `retmax`    | `int`          | Maximum records to retrieve per query         | `10000`      |
| `db`        | `str`          | NCBI database to search                        | `'nucleotide'`|
| `fmt`       | `str`          | Download format                                | `'fasta'`    |

#### Methods

##### `search(query: str) → List[str]`

Search NCBI database and return list of sequence IDs.

```python
>>> db = Database(topics=[], out_dir="./data", email="user@email.com")
>>> ids = db.search("BRCA1[Gene] AND Homo sapiens[Organism]")
>>> len(ids)
1250
```

##### `build(with_index: bool = True)`

Execute the complete download pipeline for all topics.

```python
>>> db.build(with_index=True)
[1/3] Processing topic: BRCA1[Gene] AND Homo sapiens[Organism]
Found 1250 sequences for query: BRCA1[Gene] AND Homo sapiens[Organism]
Fetching batch 1/3
Wrote 1250 sequences to ./data/BRCA1_Gene_AND_Homo_sapiens_Organism.fasta
Building index: ./data/BRCA1_Gene_AND_Homo_sapiens_Organism.idx
Successfully processed 3/3 topics
```

### Utility Functions

#### `create_index(input_dir: str, index_path: str = "combined.idx")`

Create a unified SQLite index from all FASTA files in a directory.

```python
from biosaic.database import create_index

# Index all FASTA files in directory
create_index("./fasta_files/", "combined_sequences.idx")
```

**Features:**
- Automatic duplicate sequence ID handling
- Progress logging for large datasets
- Robust error handling for corrupted files

#### `convert_fasta(input_dir: str, output_dir: str, mode: str = 'csv')`

Convert FASTA files to structured formats (CSV or Parquet).

```python
from biosaic.database import convert_fasta

# Convert all FASTA files to CSV
convert_fasta("./fasta_data/", "./csv_output/", mode='csv')

# Convert to Parquet format
convert_fasta("./fasta_data/", "./parquet_output/", mode='parquet')
```

**Output Schema:**
- `id`: Sequence identifier
- `name`: Full sequence description
- `length`: Original sequence length
- `sequence`: Complete sequence string

## Usage Examples

### Basic Usage

```python
from biosaic.database import Database

# Initialize database builder
db = Database(
    topics=[
        "BRCA1[Gene] AND Homo sapiens[Organism]",
        "TP53[Gene] AND Homo sapiens[Organism]",
        "CFTR[Gene] AND Homo sapiens[Organism]"
    ],
    out_dir="./genome_data/",
    email="researcher@university.edu",
    api_key="your_ncbi_api_key_here",
    max_rate=10.0,  # With API key, can use higher rates
    batch_size=1000,
    retmax=5000
)

# Download and process all topics
db.build(with_index=True)
```

### Using Predefined Queries

```python
from biosaic.database import Database
from biosaic.queries import EntrezQueries

# Load predefined queries
queries = EntrezQueries()
print(f"Available queries: {len(queries.queries)}")

# Use subset of queries
selected_queries = queries.queries[:10]  # First 10 queries

db = Database(
    topics=selected_queries,
    out_dir="./curated_genomes/",
    email="researcher@university.edu",
    retmax=2000
)

db.build()
```

### Advanced Pipeline with Conversion

```python
from biosaic.database import Database, convert_fasta, create_index

# Step 1: Download sequences
db = Database(
    topics=["COVID-19[All Fields]", "SARS-CoV-2[All Fields]"],
    out_dir="./covid_data/",
    email="researcher@university.edu",
    db="nucleotide",
    retmax=10000
)
db.build(with_index=True)

# Step 2: Create unified index
create_index("./covid_data/", "covid_combined.idx")

# Step 3: Convert to structured format
convert_fasta("./covid_data/", "./covid_csv/", mode='csv')
convert_fasta("./covid_data/", "./covid_parquet/", mode='parquet')
```

### Working with Different Databases

```python
# Protein sequences
protein_db = Database(
    topics=["insulin[Protein Name]", "hemoglobin[Protein Name]"],
    out_dir="./proteins/",
    email="researcher@university.edu",
    db="protein",  # Search protein database
    fmt="fasta"
)
protein_db.build()

# PubMed abstracts
pubmed_db = Database(
    topics=["CRISPR gene editing[Title/Abstract]"],
    out_dir="./abstracts/",
    email="researcher@university.edu",
    db="pubmed",
    fmt="abstract"
)
pubmed_db.build()
```

## Error Handling and Logging

The Database class includes comprehensive error handling:

```python
import logging

# Configure logging level
logging.basicConfig(level=logging.DEBUG)

db = Database(
    topics=["invalid_query[NonexistentField]"],
    out_dir="./test/",
    email="test@email.com"
)

try:
    db.build()
except Exception as e:
    print(f"Build failed: {e}")
    # Logs will show detailed error information
```

**Common Error Types:**
- `ValueError`: Invalid parameters or empty results
- `FileNotFoundError`: Invalid file paths
- `HTTPError`: Network or NCBI API issues
- `IncompleteRead`: Partial download recovery

## Performance Optimization

### Rate Limiting

```python
# Without API key (max 3 requests/second)
db_basic = Database(
    topics=queries,
    out_dir="./data/",
    email="user@email.com",
    max_rate=3.0
)

# With API key (up to 10 requests/second)
db_enhanced = Database(
    topics=queries,
    out_dir="./data/",
    email="user@email.com",
    api_key="your_api_key",
    max_rate=10.0
)
```

### Memory Management

```python
# Large dataset handling
large_db = Database(
    topics=["16S rRNA[All Fields]"],
    out_dir="./large_dataset/",
    email="researcher@university.edu",
    batch_size=100,  # Smaller batches for memory efficiency
    retmax=50000     # Large dataset
)
```

## File Organization

The Database class creates organized output structure:

```
output_directory/
├── topic1_sanitized_name.fasta
├── topic1_sanitized_name.idx
├── topic2_sanitized_name.fasta
├── topic2_sanitized_name.idx
└── ...
```

**Filename Sanitization:**
- Removes invalid characters
- Replaces spaces with underscores
- Ensures cross-platform compatibility

## Integration with Dataset Class

```python
from biosaic.database import Database
from biosaic.dataset import Dataset

# Step 1: Download sequences
db = Database(
    topics=["BRCA1[Gene] AND Homo sapiens[Organism]"],
    out_dir="./brca_data/",
    email="researcher@university.edu"
)
db.build()

# Step 2: Create dataset for training
dataset = Dataset(
    mode="dna",
    kmer=3,
    index_path="./brca_data/BRCA1_Gene_AND_Homo_sapiens_Organism.idx",
    continuous=True
)

# Step 3: Generate training batches
batch = dataset.get_batch("train", batch_size=32, block_size=512)
```

## Requirements

* **Python**: 3.8+
* **Biopython**: ≥1.78
* **pandas**: ≥1.3.0
* **pyarrow**: ≥5.0.0 (for Parquet support)
* **sqlite3**: Built-in with Python
* **NCBI Account**: Required for API key (recommended)

## Best Practices

1. **Always provide an email address** - Required by NCBI
2. **Use API keys** - Enables higher rate limits and better reliability
3. **Monitor rate limits** - Respect NCBI's usage policies
4. **Handle large datasets** - Use appropriate batch sizes
5. **Validate queries** - Test queries on small datasets first
6. **Backup data** - Keep copies of downloaded sequences

## Troubleshooting

### Common Issues

**Empty Results:**
```python
# Check query syntax
ids = db.search("your_query_here")
if not ids:
    print("No sequences found. Check query syntax.")
```

**Rate Limit Errors:**
```python
# Reduce rate and increase batch delays
db = Database(
    topics=queries,
    out_dir="./data/",
    email="user@email.com",
    max_rate=1.0  # Very conservative
)
```

**Memory Issues:**
```python
# Use smaller batch sizes
db = Database(
    topics=queries,
    out_dir="./data/",
    email="user@email.com",
    batch_size=50  # Smaller batches
)
```
