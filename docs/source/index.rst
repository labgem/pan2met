.. pan2met documentation master file, created by
   sphinx-quickstart on Mon Dec  8 14:36:19 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pan2met: Predicting metabolic networks of procaryotes at pangenome scale
=======================================================================

**pan2met** is a software suite aiming to predict metabolic network of procaryotes, bacteria and archaea, at the pangenome scale.

**pan2met** requires a pangenome in `PPanGGOLiN <https://github.com/labgem/PPanGGOLiN>`_'s HDF5 format.
It starts to annotate pangenome gene families reference proteins with enzymatic activity by homology with proteins from reference databases (KEGG, MetaCyc, UniProt). Then, It tries to predict what metabolic pathways constitute the target organism metabolism, both at the pangenome scale and at the scale of the strain.



User guide
----------

.. toctree::

   Tutorials <tutorial/index>

Developper guide
----------------

.. toctree::

   API <pydoc/modules>
