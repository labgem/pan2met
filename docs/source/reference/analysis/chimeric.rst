Chimeric pathways
=================


**pan2met** allows to perform metabolic pathway predictions for a set of genes of a pangenome.
It may happen that those prediction contains *chimeric* pathways. A chimeric pathway is a pathway composed of reactions whose enzyme are found in distinct strain in the pangenome. For instance, we could have a chimeric pathway with reactions 1, 2 and 3 with reaction 1 in a strain A and reactions 2 and 3 only in strain B.

Such a prediction is kind of an artefact, which might both reflect a lack of annotation in both strain A and B, but could also simply be an erroneous prediction where neither strain A nor strain B have the pathway.

To identify such spurious cases, **pan2met** tool suite includes a script ``chimeric.py``.
To use this script, you can use the following command:

.. code:: bash

    python3 -m pan2met.analysis.chimeric \
          --gene-presence-absence gene_presence_absence.Rtab \
          --gene-reaction gene_reactions.tsv \
          --output ./chimeric_pathways.list


``--gene-presence-absence`` expects a binary matrix in scoary .Rtab format, with genes in rows, strain in columns and 1 if the gene is present in the strain, 0 otherwise.


``--gene-reaction`` expects a TSV file with two columns. The first columns is the name of the gene family. The second column is the name of the reaction the gene family catalyzes.

The output file ``chimeric_pathways.list`` is a list of chimeric pathways.

