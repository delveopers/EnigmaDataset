# EnigmaDB

Dataset generation pipeline for Enigma2 using NCBI database. Along with additional helper such as ``Dataset`` class to retireve & create batches for training ML models.

## Prerequisites

Before setting up EnigmaDB, ensure that you have the following prerequisites installed:

1. Python 3.8 or higher
2. pip (Python package installer)

### Dependencies

- [Biopython](https://biopython.org/)
- [NumPy](https://numpy.org/)
- [PyTorch](https://pytorch.org/)
- [Biosaic](https://pypi.org/project/biosaic/)

## Installation

#### 1. From PyPI

  ```cmd
  pip install enigmadatabase
  ```

#### 2. Clone the Repo

  ```bash
  git clone https://github.com/delveopers/EnigmaDataset.git
  cd EnigmaDB
  ```

## Documentation & Usage

### Data gathering pipeline

```python
from EnigmaDB import Database, EntrezQueries
queries = EntrezQueries()   # get queries

db = Database(topics=queries(), out_dir="./data/raw", email=EMAIL, api_key=API_KEY, retmax=1500, max_rate=10)   # set parameters
db.build(with_index=False)  # startbuilding
```

### Creating Indexes

```python
from EnigmaDB import create_index

create_index("./data/raw")    # add path to data
```

### Converting versions

```python
from EnigmaDB import convert_fasta

convert_fasta(input_dir="./data/raw", output_dir="./data/parquet", mode='parquet')  # for parquet
convert_fasta(input_dir="./data/raw", output_dir="./data/parquet", mode='csv')  # for csv
```

for more technical information, refer to documentation:

1. [Database.md](https://github.com/delveopers/EnigmaDataset/blob/main/docs/Database.md)
2. [Dataset.md](https://github.com/delveopers/EnigmaDataset/blob/main/docs/Dataset.md)

## Project Structure

```text
├── docs/
├── ├── Database.md
├── ├── Dataset.md
├── src/
├── ├── __init__.py
├── ├── _database.py    # ``Database`` class for downloading data from NCBI
├── ├── _dataset.py     # ``Dataset`` a dataloader class for enigma2
├── ├── _queries.py     # contains queries for DB pipeline
├── README.md
├── setup.py
├── pyproject.toml
├── requirements.txt  # List of Python dependencies
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

1. Fork the repository.

2. Create a feature branch:

  ```bash
  git checkout -b feature-name
  ```

3. Commit your changes:
  
  ```bash
  git commit -m "Add feature"
  ```

4. Push to the branch:

  ```bash
   git push origin feature-name
  ```

5. Create a pull request.

Please make sure to update tests as appropriate.

## License

MIT License. Check out [License](LICENSE) for more info.
