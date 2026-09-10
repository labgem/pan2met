# `pan2met`

> A python library / CLI to predict metabolic networks at the pangenome scale.

## Usage

To predict metabolic pathways with `pan2met`, you will need a set of catalyzed reactions.
To identify such a set of reactions, you can use the [nextflow](https://nextflow.io) workflow [pan2met-wf](https://github.com/labgem/pan2met-wf/). This workflow supports multiple enzyme annotation sources to map proteins to MetaCyc reactions.

Given a file `reaction.list` with a list of MetaCyc reaction identifiers, you can use the following command to predict pathways using a set of decision rules.

``` bash
python3 -m pan2met metabolism \
    --reactions reaction.list \
    --output pathway.list \
    --taxon-id 561 \
    --reason pan2met.log
```


```text
usage: pan2met metabolism [-h] -r REACTIONS -o OUTPUT [--reason REASON] [-t TAXON_ID]

infer the (pan)metabolism (i.e., a set of expected metabolic pathways)

options:
  -h, --help            show this help message and exit
  -r, --reactions REACTIONS
                        a file listing the reactions found in the (pan)-reactome
  -o, --output OUTPUT   the path of the output file listing all identifiers of pathway infered to be present
  --reason REASON       the path to an output file with a reason log.
  -t, --taxon-id TAXON_ID
                        the NCBI-Taxonomy tax id of the target organism.
```

For more information on how to use `pan2met`, please refer to the [`pan2met` documentation](https://pan2met.readthedocs.io/en/latest/index.html)

## Installation

### From source

1. Clone this repository

   ``` bash
   git clone https://github.com/labgem/pan2met.git
   cd pan2met
   ```
2. Install locally in a pixi virtual environment

   `pan2met` uses [graph-tool](https://graph-tool.skewed.de) to manage the pangenome graph datastructure.
   As `graph-tool` is not available on [PyPI.org](https://pypi.org/), being a C++ backed Python package, you will need to install graph-tool on your own with your OS package management system, or using the conda-forge distribution.
   For an easy environment creation, you can use [pixi](https://prefix.dev/):

   ``` bash
   pixi shell
   ```

   Then, the `pan2met` command line interface should be installed and available in your PATH:

   ``` bash
   pan2met --version
   ```

## Setup

Create a configuration file, in e.g. `conf/configuration.ini`, from [provided example configuration file](./src/pan2met/conf/default.ini).

You will most probably need to adapt the `[reference]` section.
1. Update `ncbi_taxonomy` directory path, with the path to the directory with a dump of the NCBI-Taxonomy.
2. Update the `source` key, to either `metabiantes` or `padmet` depending of the format of metabolism knowledge base to use.

If you use `metabiantes` as the reference knowledge base for metabolism, please refer to the [metabiantes git repository](https://github.com/labgem/metabiantes/) for instructions on how to setup a `metabiantes` SQL database.

### How to download support materials

#### Reference taxonomy from the NCBI Taxonomy database

Download a NCBI Taxonomy dump from <https://ftp.ncbi.nih.gov/pub/taxonomy/> to a local folder, and adapt the path in the configuration file in section `[reference]`, key `ncbi_taxonomy`.
