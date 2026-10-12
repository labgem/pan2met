Per genome pathway completion
=============================

As **pan2met** can also be used on a genome, instead of pangenome, we added a script to compute the completion value of metabolic pathways for a genome.

In this mode, we need a folder with the sets of reactions that may be catalyzed in the cells of the target organism.

We use the basename of the file as the identifier of the target genome in the output completion matrix table.

You can use this script as follows:

.. code:: bash

    python3 -m pan2met.analysis.completion_genomes \
          --genome-reactions ./folder/ \
          --pathways ./pathways.list \
          --prefix "completion_"

Contrary to the pangenome case (:doc:`completion_pangenome`), we do not have information on the pangenome partition, nor global completion values.

The output completion table columns are:

.. csv-table:: Columns of the pathway completion matrix
   :header: "Column", "Content"
   :widths: 15, 10


            pathway,"the identifier of the pathway"
            pathway name,"the name of the pathway, if available in the knowledgebase"
            nb reactions,"number of reactions in the pathway"
            "max completion","maximum strain-wise completion"
            "genome completion","completion of the pathway for the given genome"
