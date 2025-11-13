# Pangenome to Panmetabolome

A python library / CLI to infer metabolic networks at the pangenome scale.

## Installation

### From source

1. Clone this repository

``` bash
git clone https://gitlab.com/sortion/pangenome2panmetabolome.git
cd pangenome2panmetabolome
```

2. Install locally in a virtual environment

``` bash
python3 -m venv .venv/pangenome2panmetabolome
source .venv/pangenome2panmetabolome/bin/activate
pip install -e .
```

## Configure database backend

#### SQLite3 backend for BioPAX RDFs

1. Install [rdftab.rs](https://github.com/ontodev/rdftab.rs).
2. Export MetaCyc BioPAX level 3 OWL file.
3. Replace biocyc.org with biopax.org in the owl file:
   ``` bash
   sed 's$http://biocyc.org/biopax/biopax-level3\#$http://www.biopax.org/release/biopax-level3.owl\#$g' metacyc-biopax-level3.owl > metacyc-biopax-level3-biopax-prefix.owl 
   ```
4. Prepare a SQLite database with a table of RDF prefix:
   ``` bash
   sqlite3 metacyc-biopax.db < resources/prefix.sql
   ```
5. Load the RDF from the OWL file into the SQLite database:
   ``` bash
   rdftab metacyc-biopax.db < metacyc-biopax-level3.owl
   ```
