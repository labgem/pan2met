# `pan2met`

> A python library / CLI to infer metabolic networks at the pangenome scale.

## Installation

### From source

1. Clone this repository

   ``` bash
   git clone https://gitlab.com/sortion/pan2met.git
   cd pan2met
   ```
2. Install locally in a virtual environment

   ``` bash
   python3 -m venv .venv/pan2met
   source .venv/pan2met/bin/activate
   pip install -e .
   ```

## Setup

Create a configuration file, in e.g. `conf/configuration.ini`, from provided example configuration file:

``` ini
[reactome]

[inference]
pathway_score_cutoff = 0.35

[reference]
ncbi_taxonomy = /mnt/shared/bank/NCBI-Taxonomy/taxdmp_2026-01-01
source = metabiantes

[metabiantes]
database = metabiantes
```
You will most probably need to adapt the `[reference]` section.
1. Update `ncbi_taxonomy` directory path, with the directory where lies the NCBI-Taxonomy dump.
2. Update the `source` key, to either `metabiantes` or `padmet` depending of the _source_ of metabolism knowledge.

If you use `metabiantes` as the reference knowledge base for metabolism, please refer to [metabiantes git repository](https://gitlab.com/sortion/metabiantes/) for instructions on how to setup a `metabiantes` SQL database.

### Download support materials

#### Reference taxonomy from the NCBI-Taxonomy database

Download a NCBI Taxonomy dump from <https://ftp.ncbi.nih.gov/pub/taxonomy/> to a local folder, and adapt the path in the configuration file in section `[reference]`, key `ncbi_taxonomy`.

## Fix clyngor clingo Answer set parsing error

On clingo v5.8.0 (at least), the output Answer set is followed by the execution time, e.g.:

``` text
Answer 1 (Time: 0.608s)
```

This cause an int parse error in clyngor (see [clyngor merged pr #34](https://github.com/Aluriak/clyngor/pull/34)), until a new version is available on PyPI including this fix, if you have a clingo version having the execution time in answer set output, you should install clyngor with:

``` bash
pip install git+https://github.com/Aluriak/clyngor@master
```
