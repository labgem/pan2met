How to infer a "minimal" set of metabolic pathway with a MinPath-like approach
==============================================================================


The MinPath algorithm, from Ye and Doak 2009 paper "A Parsimony Approach to Biological Pathway Reconstruction/Inference for Genomes and Metagenomes" (`doi <https://doi.org/10.1371/journal.pcbi.1000465>`__) relies on a Mixed Integer Linear programming approach to select a minimum set of metabolic pathway that covers a set of target reactions.

This approach can be reduced to finding a solution to a `set cover <https://en.wikipedia.org/wiki/Set_cover_problem>` problem, where the universe is the whole set of reactions in the knowledgebase, the target set we want to cover are the reaction catalyzis found in the genome and candidate subsets are the sets of reactions for each metabolic pathways.

In term of Mixed Integer Linear programming, as proposed by Ye and Doak, this problem can be expressed as:

.. math::

   \min \sum_{j=1}^p P_j

   s.t. \sum_{j=1}^p M_{ij} P_j \geq 1 & \forall i \in [1, n]


where :math:`P_j` equals 1 when the pathway :math:`j` is selected, 0 otherwise; and :math:`M_{ij}` equals 1 when the pathway :math:`j` contains the reaction :math:`i`.

We can also solve the set cover problem with Answer Set Programming, with the following clingo sketch (src/pan2met/asp/rules/minpath_like.lp):

.. literalinclude:: ../../../src/pan2met/asp/rules/minpath_like.lp


To use this ASP model to predict a minimal subset of pathway covering a set of reactions, you can use **pan2met**'s subcommand `minpath`:

.. code:: bash

    pan2met minpath --reactions reaction.list --output pathways.list

Beware that if a reaction *cannot* be covered by any pathway (when it is not represented in any pathway in the knowledgebase, for instance), the ASP model will be unsatisfiable and no minimal list of pathway can be returned.
