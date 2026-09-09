# `pan2met`

> A python library / CLI to predict metabolic networks at the pangenome scale.

## Usage

To predict metabolic pathways with `pan2met`, you will need a set of catalyzed reactions.
To identify such a set of reactions, you can use the [nextflow](https://nextflow.io) workflow [pan2met-wf](https://github.com/labgem/pan2met-wf/). This workflow supports multiple enzyme annotation sources to map proteins to MetaCyc reactions.

Given a file `reaction.list` with a list of MetaCyc reaction identifiers, you can use the following command to predict pathways using a set of decision rules.

``` bash
python3 -m pan2met.inference.pathologic.pythonic \
    --reactome reaction.list \
    --output pathway.list \
    --taxon 561 \
    --reason pan2met.log
```

The `pan2met.inference.pathologic.pythonic` module command line interface requires the following options:

- `--reactome` -- a file with a list of reaction identifiers
- `--output` -- output file with the list of predicted pathways
- `--taxon` -- an integer corresponding to the taxonomy identifier from the NCBI Taxonomy
- `--reason` -- a file to log a reason leading to keep the pathway or to reject it


## Installation

### From source

1. Clone this repository

   ``` bash
   git clone https://github.com/labgem/pan2met.git
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
pathway_score_threshold = 0.35

[reference]
ncbi_taxonomy = /mnt/shared/bank/NCBI-Taxonomy/taxdmp_2026-01-01
knowledge_base = metabiantes

[metabiantes]
database = /mnt/shared/bank/
```
You will most probably need to adapt the `[reference]` section.
1. Update `ncbi_taxonomy` directory path, with the path to the directory with a dump of the NCBI-Taxonomy.
2. Update the `source` key, to either `metabiantes` or `padmet` depending of the format of metabolism knowledge base to use.

If you use `metabiantes` as the reference knowledge base for metabolism, please refer to [metabiantes git repository](https://github.com/labgem/metabiantes/) for instructions on how to setup a `metabiantes` SQL database.

### Download support materials

#### Reference taxonomy from the NCBI Taxonomy database

Download a NCBI Taxonomy dump from <https://ftp.ncbi.nih.gov/pub/taxonomy/> to a local folder, and adapt the path in the configuration file in section `[reference]`, key `ncbi_taxonomy`.
