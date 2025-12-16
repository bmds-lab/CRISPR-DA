# CRISPR-DA
Adapting CRISPR gRNA design for detection assays.

## Preamble

>CRISPR-Cas systems offer a viable alternative to traditional detection and diagnosis methods. However, their effectiveness relies heavily on the selection of appropriate guide RNA sequences. Existing gRNA design tools were primarily developed for gene editing and are not always directly applicable to CRISPR-based detection assays. In particular, alignment-based methods are still used to estimate gRNA specificity, even though they can miss a substantial portion of off-target sites. In this work, we introduce CRISPR-DA, a CRISPR gRNA design tool for detection assays. We show that it provides a better assessment of gRNA specificity than BLAST, which detected only 33.27 % and 0.43 % of cross-species off-targets in two datasets. Additionally, CRISPR-DA ran two and six times faster than BLAST on these datasets, respectively. Our method incorporates advances from gene-editing guide RNA design tools, including uncertainty-informed guide RNA design, to improve the selection of guides with high on-target activity. CRISPR-DA is available at https://github.com/bmds-lab/CRISPR-DA.

## Licence
BSD 3-Clause "New" or "Revised" License

See `LICENCE` and https://opensource.org/license/BSD-3-Clause for more details.

## Dependencies

### Option 1. Using Anaconda (Recommended)
1. Download the `environment.yml` file from GitHub.

2. Download and install [Anaconda](https://www.anaconda.com/download).

3. Setup conda enviroment.
   ```
   conda env create -f <path/to/downloaded/environment.yml>
   ```

### Option 2. Install dependencies manually
1. Download and install the following (The software/program version provided were using during testing but there "should" be no problem using newer version).
-  [Python](https://www.python.org/)>=3.12
-  [CMake](https://cmake.org/)>=4.1.2
-  [Boost](https://www.boost.org/)>=1.82
-  [Make](https://www.gnu.org/software/make/)>=4.2.1
NOTE: You can replace Make with whichever build system you like (E.g. [Ninja](https://ninja-build.org/)).

## Installation

### From [PyPi](https://pypi.org/project/crispr-da/) using pip (Recommended)
1. Install using pip.
   ```
   pip install crispr-da
   ```

2. Run `config` to install required binaries and initialise configuration file.
   ```
   crispr_da config
   ```

### From source
1. Clone the repository from GitHub.
   ```
   git clone https://github.com/bmds-lab/CRISPR-DA
   ```

2. Build distribution files.
   ```
   python3 -m build
   ```
   NOTE: It may be useful to ensure the latest version of PyPA's build is install. Simply run `python3 -m pip install --upgrade build`.

3. Install CRISPR-DA using distribution files and pip.
   ```
   python3 -m pip install ./dist/*.whl
   ```

4. Run `config` to install required binaries and initialise configuration file.
   ```
   crispr_da config
   ```

## CRISPR-DA CLI

### Usage

   ```
   $ crispr_da -h
   usage: crispr_da [-h] [-v] {config,analyse} ...

   CRISPR-DA: gRNA design for detection assays.

   options:
   -h, --help        show this help message and exit
   -v, --version     Print CRISPR-DA version

   subcommands:
   {config,analyse}
      config          Run config
      analyse         Run analysis
   ```

### Example

   ```
   $ crispr_da analyse --target_gene_id 43740568 --evaluation_accessions GCA_000820495.2 GCA_000838265.1
   ```

Expected output format (minimised for this readme):

| accession       | TACACTAATTCTTTCACACGTGG |  ...  | CTTTCTTTTCCAATGTTACTTGG |
|:----------------|:------------------------|:------|:------------------------|
| GCA_000820495.2 | (100.0, 100.0)          |  ...  | (100.0, 100.0)          |
| GCA_000838265.1 | (100.0, 100.0)          |  ...  | (100.0, 100.0)          |

## CRISPR-DA Python Library

### Usage
CRISPR-DA can also be access using Python.
```
from crispr_da import run_analysis
```

### Example
An minimal example has been included in the `examples/design_for_covid_spike_protein.py`