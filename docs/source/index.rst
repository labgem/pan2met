.. pan2met documentation master file, created by
   sphinx-quickstart on Mon Dec  8 14:36:19 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pan2met: Predict metabolic networks of procaryotes at pangenome scale
========================================================================

**pan2met** is a software suite aiming to provide tools to reconstruct metabolic networks of procaryotes: bacteria and archaea, at the pangenome scale.

**pan2met** can use a pangenome gene graph, such as the ones produced by `PPanGGOLiN <https://github.com/labgem/PPanGGOLiN>`_ or `panaroo <https://github.com/gtonkinhill/panaroo>`_ to try to enhance the metabolic pathway predictions.

The decision rules implemented in **pan2met** can also be used to reconstruct a metabolic network from a single genome.

In spite of being designed primary with prokaryote genomes and pangenomes in mind, **pan2met** could also be applied to a certain extent to Eukaryotes.


User guide
----------

.. toctree::

   Installation <setup>
   Tutorials <tutorial/index>
   Reference <reference/index>


Developer guide
----------------

.. toctree::

   API <pydoc/modules>
