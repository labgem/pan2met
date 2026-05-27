Analyses on the results of pan2met
==================================

Completion of the metabolic pathways per strain
-----------------------------------------------


Identification of chimeric pathways
-----------------------------------

If pan2met is used on a set of enzymes identified from a pangenome, the reactions of the predicted metabolic pathways may be catalyzed by enzymes of certain strain of the pangenome. This leads to what we call _chimeric_ pathways.

The Python script analysis.chimeric enable to identify such pathways.

This script takes a .Rtab file from PPanGGOLiN with a gene presence absence matrix, a set of predicted pathways and a mapping from gene/protein families to reaction identifiers and returns the list of chimeric pathways. Pathways are considered chimeric if any pair of reactions is catalyzed in a completely distinct pair of set of strains, for reaction catalyzed in at least one strain.

An example of command line to run such an analysis is the following:

.. code:: bash

    python3 -m pan2met.analysis.chimeric \
    --pathways pathways.list \
    --gene-presence-absence ppanggolin/gene_presence_absence.Rtab \
    --gene-reaction gene_reaction.tsv


The file pathways.list is a list of pathway identifiers, one pathway identifier per row.

The file gene_presence_absence.Rtab is a binary matrix in .Rtab format, with gene family identifiers as row index and strain identifiers as column index. A value as coordinate $(i, j)$ is a Boolean value, 1 if the gene family $i$ is found in strain $j$.

The file gene_reaction.tsv is a two column TSV file, with column 1 corresponding to gene identifiers, column 2 to reaction identifiers, for reaction catalyzed by the enzyme coded by the gene.

The analysis is done on the set of non-orphan, non-spontaneous reactions of the pathway.

You can also provide a config file with option --config to override the default config, following the syntax of :doc:`config`.
