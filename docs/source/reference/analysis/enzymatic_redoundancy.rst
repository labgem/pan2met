Enzymatic redoundancy
=====================

A reaction can be catalyzed by several enzymes in a strain, or by different enzymes in different strains that composes the pangenome.

**pan2met** includes a script to compute the number of enzyme catalyzing a reaction.

We named this script ``enzymatic_redoundancy.py``.

You can use this script as follows:


.. code:: bash

    python3 -m pan2met.analysis.enzymatic_redoundancy \
        --pangenome pangenome.h5 \
        --gene-reaction gene_reaction.tsv \
        --output reaction_enzyme_count.tsv

The output table ``reaction_enzyme_count.tsv`` contains the following two columns:
- reaction
- enzyme count


This script supports also the option ``--by-strain`` to count the number of enzymes catalyzing the reaction, within each strain.
When this flag is in use, the output table columns becomes:
- reaction
- strain
- enzyme count

