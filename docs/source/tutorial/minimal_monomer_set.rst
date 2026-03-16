Infer the minimal set of monomers required to catalyzes a set of reactions
==========================================================================

`pan2met` includes a method to compute an "inverse" problem of the reactome.

We define the "inverse" problem of the reactome as follows:

:Given: A set of reactions and the catalysis rules

:Return: A (minimal) set of monomer that enable the catalyzes of all these reactions.



The method implemented in `pan2met` uses Clingo Answer Set Programming solver.

It proceeds in two steps:

1. Find the set of monomer that may be involved in the catalysis of the reactions.
2. Find the minimal subset of such a set of monomer.


Currently, it is only possible to use the MetaCyc knowledge base as a source of Gene - Protein - Reaction rules, but it should be possible to use other sources of data.

Preparation step
----------------

The first step to do this analysis is to construct two files:
1. A 'forward' GPR rule, containing rules such as:

   .. code:: prolog

      complex("CPLX-1") :- monomer("MONOMER-A") , monomer("MONOMER-B").

   meaning CPLX-1 complex is composed of both MONOMER-A and MONOMER-B,

   .. code:: prolog

      reaction("RXN-1") :- monomer("MONOMER-C").
      reaction("RXN-2") :- complex("CPLX-1").

    meaning reaction RXN-1 is catalyzed by monomer MONOMER-C and inferred to be in the reactome whenever MONOMER-C is in the proteome, and similarly for reaction RXN-2.

2. A reverse GPR rule, similar to the forward one, except `complex` atoms are renamed `potential_complex`, `monomer` renamed `potential_monomer`, and the head and tail of ASP rules are reversed, and whenever the tail contains a conjonction of atoms, each of these atom will be a head of one of the converted rule.


To generate these files:

First, launch pathway-tools python API:

.. code:: bash

  pathway-tools -lisp -python-local-only-non-strict

Then, on another shell:

.. code:: bash

  python3 -m src.pan2met.asp.gpr_rules -o ./tmp/metacyc_gpr.lp

  python3 -m src.pan2met.asp.reverse_gpr_rules -i ./tmp/metacyc_gpr_rules -o ./tmp/metacyc_reverse_gpr.lp

Application on an example
-------------------------

In `resources/tests/reactome_inference`, there are sample .lp files to test the approach.

First, generate the list of potential monomer, namely, monomer that can contribute to the catalysis of the reactions:

Construct a file with some atoms `reactions/1` `reactions.lp`:

.. literalinclude:: ../../../resources/tests/reactome_inference/reactions.lp

Call clingo to infer the set of potential monomers:

.. code:: bash

  clingo -V0 --out-atomf=%s. ./tmp/metacyc_reverse_gpr.lp ./resources/tests/reactome_inference/reactions.lp <(echo "#show potential_monomer/1.") | head -1 > /tmp/potential_monomers.lp

The command `clingo -V0 --out-atomf=%s.` ensures that only the atoms are outputed. The `<(echo "#show potential_monomer/1.")` file input insert a `#show` directive to keep potential_monomer only, and not the `reaction/1` atoms, if they would appear in the answer set, the next step would not be working, as `reaction` would already be satisfied.

The output file will look like:

.. literalinclude:: ../../../resources/tests/reactome_inference/potential_monomers.lp

Then, construct a file `target_reactions.lp`, with the same reactions identifiers as in `reactions.lp`, except the predicate `target_reaction` is used instead of `reaction`:

.. literalinclude:: ../../../resources/tests/reactome_inference/target_reactions.lp

Finally, run clingo once more:

.. code:: bash

  clingo ./tmp/metacyc_gpr.lp ./src/asp/required_monomer_given_reactions.lp /tmp/potential_monomers.lp ./resources/tests/reactome_inference/target_reactions.lp


The ASP program `./src/asp/required_monomer_give_reactions.lp` contains the logic to choose a subset of the potential_monomers that satisfies the target set of reactions.

The selected monomers are represented by the predicates `selected_monomer/1`.
