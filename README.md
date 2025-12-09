# CRISPR-DA
CRISPR-Cas9 sgRNA design across many genomes

## Installation

1. Clone the repository from GitHub

   ```
   git clone https://github.com/bmds-lab/CRISPR-DA
   ```
   
2. Install using pip

   ```
   cd CRISPR-DA
   python3.8 -m pip install -e .
   ```
   
3. Run an example

    ```
    cd examples
    python3.8 design_for_covid_spike_protein.py
    ```

## Utilities

### CRISPR-DA CLI

#### Usage

   ```
$ CRISPR_DA --help
usage: CRISPR_DA [-h] -g GENEID -a ACCESSIONS [ACCESSIONS ...] [-o OUTPUT]

CRISPR-DA: for designing CRISPR-Cas9 sgRNA when many genomes are to be considered.

optional arguments:
  -h, --help            show this help message and exit
  -g GENEID, --geneid GENEID
                        The NCBI gene ID to extract potential CRISPR target sites from
  -a ACCESSIONS [ACCESSIONS ...], --accessions ACCESSIONS [ACCESSIONS ...]
                        A list of NCBI accessions to score CRISPR sites against
  -o OUTPUT, --output OUTPUT
                        Output files prefix (e.g. `/home/user/results/covid` would become `/home/user/results/covid-guides.csv`, etc.). Default is
                        `--geneid`.
```

#### Example

   ```
   $ CRISPR_DA -g 43740568 -a GCA_000820495.2 GCA_000838265.1 -o spike
   ```

Output file `spike-guides.md` (minimised for this readme):


| seq                     |   start |   end | ssline_1                   | ssStructure                | strand   |   acceptedByMm10db |   passedG20 |   acceptedBySgRnaScorer |   consensusCount |   ssEnergy |   passedAvoidLeadingT | passedSecondaryStructure   |   AT |   sgrnascorer2score |   passedATPercent |   passedTTTT |
|:------------------------|--------:|------:|:---------------------------|:---------------------------|:---------|-------------------:|------------:|------------------------:|-----------------:|-----------:|----------------------:|:---------------------------|-----:|--------------------:|------------------:|-------------:|
| TACACTAATTCTTTCACACGTGG |      81 |   104 | GACACUAAUUCUUUCACAGUGCUUUU | ..((((..(((......(((().... | +        |                  0 |           1 |                       1 |                2 |      -19.2 |                     0 | !                          |   65 |           1.00437   |                 1 |            1 |
| CTCAGTTTTACATTCAACTCAGG |     134 |   157 | GUCAGUUUUACAUUCAACGUGCUUUU | ..((.((((........(((.).... | +        |                  0 |           0 |                       0 |                0 |      -19.2 |                     1 | 0                          |   65 |          -0.537641  |                 1 |            0 |



Output file `spike-scores.md` (minimised for this readme):


| accession       | TACACTAATTCTTTCACACG   | CTGGGACCAATGGTACTAAG   | GAGAAGTCTAACATAATAAG   | GAAGTCAGACAAATCGCTCC   |
|:----------------|:-----------------------|:-----------------------|:-----------------------|:-----------------------|
| GCA_000820495.2 | (100.0, 100.0)         | (100.0, 100.0)         | (100.0, 100.0)         | (100.0, 100.0)         |
| GCA_000838265.1 | (100.0, 100.0)         | (100.0, 100.0)         | (100.0, 100.0)         | (100.0, 100.0)         |


## Dependencies

- Python 3.7 +

- [Pandas](https://pypi.org/project/pandas/)

- [NCBI Dataset PyLib](https://pypi.org/project/ncbi-datasets-pylib/)

