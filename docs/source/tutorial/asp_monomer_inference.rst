Using an Answer-Set-Programming approach to identify pathway, from Monomers to pathways through reactions and protein complexes
===============================================================================================================================

`pan2met` implements a basic/naive approach for the pathway prediction problem, that can serve as a base for further method refinement.

The rule is simple: when a protein monomer catalysing a reaction, is found to have an homologous protein in the target organism proteome: infer the reaction to be present. When every monomer component of a protein complex is present in an cell, infer the presence of the protein complex. Further populate the set of potential enzymes and reactions with this proteins.

Finally, given the set of inferred reactions, simply infer all metabolic pathway with all reactions present in the reactome to be also present in the metabolome of the target organism.


Preparation step
----------------

A prerequisite for such an Answer-Set-Programming based inference is first to encode the knowledge base with AnsProlog atoms and rules.
To do so, `pan2met` offers a subcommand named TODO.
