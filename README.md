# Pangenome to Panmetabolome

A python library / CLI to infer metabolic networks at the pangenome scale.

## Installation

### Requirements
To install gmpy library of PPanGGOLiN, you need libgmp "GNU Multiple Precision library",  GNU MPFR library

On debian install `libgmp-dev`, ...

On fedora install `mpfr-devel`, `gmp-devel`, `libmpc-devel`, `blosc-devel`.
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

## Metabolism knowledge base configuration

There are two alternatives currently under exploration to serve as a basis for the pathway knowledge base.


### Configure a SQLite backend for BioPAX RDFs

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

### Configure a SPARQL query endpoint with BioPAX RDFs

Requirements: docker or podman, and podman / docker-compose.

1. Launch the Apache Fuseki SPARQL Query endpoint, using the provided `container/compose.yaml` file; note the credentials to connect on the web interface.
   ``` bash
   pushd container
   podman compose up
   popd
   ```
2. Import the BioPAX OWL file into the Fuseki RDF triplestore database, using the web interface.
3. Copy the .env.example file to .env and adapt the Fuseki API credentials.



## Gene-Protein-Reaction rules

`pangenome2panmetabolome` includes a script to export Gene-Protein-Reaction rules in Answer Set Programming (AnsProlog) format from MetaCyc database using pythoncyc.
To do so, make sure to install [pythoncyc](https://github.com/networkbiolab/PythonCyc).

Then, launch PathwayTools python API server:
``` bash
pathway-tools -lisp -python-local-only-non-strict
```
And, launch the export script:
``` bash
python3 -m src.pangenome2panmetabolome.gpr --output <OUTPUT>.lp
```

The scripts writes two kinds of rules:

For every monomer enzyme catalyzing a reaction:
``` prolog
reaction("RXN-ID") :- monomer("MONOMER-ID").
```
For every complex protein catalyzing a reaction:
``` prolog
reaction("RXN-ID") :- complex("CPLX-ID").
```
Following these declarations, the atoms reaction/1 will be true whenever any of the atoms of the possible enzymes that catalyzes the reaction is true. This can be translated also as a rule with "OR": if "MONOMER-ID" or "CPLX-ID" (provided MONOMER-ID and CPLX-ID corresponds to enzyme catalyzing the reaction), then the reaction is inferred.

For the protein complexes, for instance a complex CPLX-1 composed of MONOMER-A and MONOMER-B, the rule is

``` prolog
complex("CPLX-1") :- monomer("MONOMER-A") , monomer("MONOMER-B").
```
where the comma, in AnsProlog, corresponds to the logical AND.


Based on these kind of rule, It is easy to infer the presence of a reaction in a reactome based solely on the presence of monomers.

For instance:

``` prolog
% GPR rules
reaction("RXN-20780") :- complex("CPLX-9428").
complex("CPLX-9428") :- monomer("G185E-7504-MONOMER") , monomer("G185E-7503-MONOMER").

% Seed
monomer("G185E-7504-MONOMER").
monomer("G185E-7503-MONOMER").

#show reaction/1.
```

``` text
Answer: 1
reaction("RXN-20780")
SATISFIABLE
```



## Fix clyngor clingo Answer set parsing error

On clingo v5.8.0 (at least), the output Answer set is followed by the execution time, e.g.:

``` text
Answer 1 (Time: 0.608s)
```

This cause an int parse error in clyngor (see [clyngor merged pr #34](https://github.com/Aluriak/clyngor/pull/34)), until a new version is available on PyPI including this fix, if you have a clingo version having the execution time in answer set output, you should install clyngor with:

``` bash
pip install git+https://github.com/Aluriak/clyngor@master
```


## NCBI-Taxonomy

Download NCBI Taxonomy dump in https://ftp.ncbi.nih.gov/pub/taxonomy/
