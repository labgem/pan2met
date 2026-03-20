Howto: Using the PathoLogic-like decision rule algorithm: two alternative implementations
=========================================================================================

The first implementation is based on a translation of the PathoLogic decision rules in Python.
The second implementation is an hybrid appproach with some PathoLogic rules implemented in Answer Set Programming (ASP) and the heuristic decision rule based on a python implementation.

Both implementations are available in **pan2met** and can be run from the command line.

Using the Answer Set Programming implementation of the PathoLogic-like decision rules
-------------------------------------------------------------------------------------

What you will need:
- A knowledge base of metabolic pathways, reactions and enzymes. In this exemple, we use the 'metabiantes' approach, a knowledge base with MetaCyc data backed by a relational database.
- A set of reactions inferred to be present in the target organism reactome.

1. Export the metabolic knowledge base as an ASP program

Using **pan2met**'s `kb_as_asp` module, we can export the metabolic knowledge base as ASP facts.

.. code:: bash

    python3 -m pan2met.asp.kb_as_asp --kb metabiantes -o ./tmp/metabiantes_as_asp.lp

2. Infer the target organism reactome

   Produce a file reactome.txt containing the list of reaction identifiers from the target organism reactome.
   This file can be derived from a genome annotation, or by sequence homology with reference monomer proteins of MetaCyc, EcoCyc, or KEGG Orthology groups.

3. Run the **pan2met** PathoLogic-like algorithm ASP implementation

   We can now run the ASP program, using the metabolic knowledge base and the target organism reactome as input.

   .. code:: bash

     python3 -m pan2met.inference.pathologic.aspic --kb-asp ./tmp/metabiantes_as_asp.lp --reactome reactome.txt --output predicted_pathways.txt

Using the python only implementation of the PathoLogic-like decision rules algorithm
------------------------------------------------------------------------------------

What you will need: (same as for the ASP implementation)
- A knowledge base of metabolic pathways, reactions and enzymes. In this exemple, we use the 'metabiantes' approach, a knowledge base with MetaCyc data backed by a relational database
- A set of reactions inferred to be present in the target organism reactome.

This time, we do not need to export the knowledge base as an ASP program, as the python implementation will directly query the knowledge base to get the necessary information for the inference at run time.

1. Infer the target organism reactome

   Produce a file reactome.txt containing the list of reaction identifiers from the target organism reactome.
   This file can be derived from a genome annotation, or by sequence homology with reference monomer proteins of MetaCyc, EcoCyc, or KEGG Orthology groups.

2. Run the **pan2met** PathoLogic-like algorithm python implementation

   We can now run the python program,

   .. code:: bash

     python3 -m pan2met.inference.pathologic.pythonic --reactome reactome.txt --taxon 562 --output predicted_pathways.txt


In the best case scenario, both implementations should give the same results. Discrepancies between the two should be taken as bug in one of the two implementations.
