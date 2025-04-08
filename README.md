# Module Gene Overlap Analysis

A web application for analyzing gene overlaps between different modules across datasets. This tool enables researchers to identify statistically significant overlaps between gene modules from their own data and reference datasets.

## Features

- Upload your own gene module data in CSV or Excel format
- Compare against built-in reference datasets
- Perform statistical testing using Fisher's Exact Test or Hypergeometric Test
- Visualize overlaps with interactive heatmaps
- Review significant overlaps with adjustable significance thresholds
- Download lists of overlapping genes for further analysis
- Apply Benjamini-Hochberg correction for multiple testing

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/module-gene-overlap.git
   cd module-gene-overlap
   ```

2. Create a virtual environment (recommended):
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install required packages:
   ```
   pip install -r requirements.txt
   ```

4. Create the data directory:
   ```
   mkdir -p data
   ```

5. Add your reference datasets (optional):
   - Place an Excel file at `data/data.xlsx` with reference datasets in separate sheets
   - Each dataset should have at least the columns 'Gene' and 'Module'

## Usage

1. Start the application:
   ```
   python server.py
   ```

2. Open a web browser and navigate to:
   ```
   http://127.0.0.1:8050/
   ```

3. Upload your dataset:
   - Your file must contain at least the columns 'Gene' and 'Module'
   - Both Excel (.xlsx) and CSV (.csv) formats are supported
   - Use the provided example file as a reference for formatting

4. Configure analysis parameters:
   - Select a reference dataset to compare against
   - Adjust the total gene universe size (default: 20,000)
   - Set the significance threshold (default: 0.05)
   - Choose minimum overlap count (default: 5)
   - Select the statistical test method

5. Run the analysis and view results:
   - Interactive heatmap shows overlap significance (-log10 of adjusted p-value)
   - Table of significant overlaps with key statistics
   - Click on any overlap to view and download detailed gene lists

## Input File Format

Your input file should contain at least the following columns:
- `Gene`: Gene identifiers (symbols or IDs)
- `Module`: Module/cluster identifiers for each gene

Example:
```
Gene,Module
BRCA1,module1
TP53,module1
MYC,module2
PTEN,module2
```

## Statistical Methods

The application offers two statistical approaches:

1. **Fisher's Exact Test**: Tests the independence of row and column variables in a contingency table.

2. **Hypergeometric Test**: Calculates the probability of finding k or more overlapping genes between two sets.

Both methods are adjusted for multiple testing using the Benjamini-Hochberg procedure to control the false discovery rate.

## Directory Structure

```
module-gene-overlap/
├── server.py               # Main application code
├── data/                # Directory for datasets
│   ├── data.xlsx        # Built-in reference datasets
│   └── example.csv      # Example dataset for users
├── requirements.txt     # Python dependencies
└── README.md            # This documentation file
```

## Requirements

- Python 3.7+
- See requirements.txt for Python package dependencies

## License

[MIT License](LICENSE)

## Citation

If you use this tool in your research, please cite:

```
Li, Huihong, et al. "Shared Transcriptomic Signatures Reveal Synaptic Pruning as a Link Between Alzheimer’s Disease and Epilepsy." bioRxiv (2024): 2024-10.
```
