Decision rules
==============

This document aims to answer the following question:
> What are the
decision rules implemented in this version of ``pan2met``?

How to enable or disable a decision rule?
-----------------------------------------

Each decision rule has an unique identifier. To enable or disable a
decision rule in the metabolic pathway inference, edit the configuration
file section ``rules`` and add or remove the corresponding identifier.

See :doc:`config` for more information on **pan2met** configuration file.


How does the inference algorithm work?
--------------------------------------

The decision algorithm will pass in review each of the decision rules.
Each of them have three possible outcomes. Either it accepts the pathway
or rejects it, or does not return a decision otherwise. In the later
case, the algorithm pass to the next decision rule. The decision rules
are ordered by an arbitrarily defined priority order. If no decision
rule accepted nor rejected the pathway at the end of the review, the
pathway is rejected by default.

A glimpse into the decision rules proposed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pathway_ontology`` - Reject some metabolic pathway based on pathway ontology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``biosynthesis_missing_last`` - Reject if a biosynthesis pathway lacks its "last" reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: If the last reaction leading to the metabolite of
interest is missing, it might be a sign that the metabolite is not
synthesized in the organism, hence that the biosynthesis pathway is not
there.

**Caution remak**: To identify the last reaction of a pathway, we use a topological sort of the reaction graph of the pathway. This assumes that the pathway is a directed acyclic graph, which might not always be the case for biosynthesis pathway. Moreover, for branching pathways, the last reaction of the obtained topological sort might be different from the reaction leading to the metabolite of interest. Thus, use this rule with caution.

``degradation_missing_first`` - Reject if a degradation pathway lacks its "first" reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: If the first reaction in a degradation pathway, the one
that initially transformed the metabolite targetted by the degradation
pathway, is not catalyzed, it might be a sign that no such degradation
occurs in the organism.

**Caution remark**: The same caution must be taken regarding the computation of the "first" reaction with topological sorting as for the rule ``biosynthesis_missing_last``.


.. ```` - Reject the transport pathways
.. ------------------------------------

``signaling_pathways`` - Reject the signaling pathways
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``all_reactions_catalyzed`` - Accept a pathway when all its non spontaneous reactions are catalyzed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: If we found evidence for the catalysis of all reactions of a pathway in an organism, we can reasonably expect this pathway to occurs in effect in the organism.

When all reactions of a given pathway are either spontaneous, or have a known enzyme catalyzing them.

When all reactions of the pathway can be spontaneous, this rule will fire too.

``all_reactions_missing`` - Reject a pathway when none of its non spontaneous reactions is catalyzed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: When none of the reactions of the pathway have an associated enzyme, we have no reason to retained it.

``KeyReaction``
^^^^^^^^^^^^^^^

``SynthesisMissingLast`` - Reject a biosynthesis pathway missing its last reaction.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


``DegradationMissingFirst`` - Reject a degradation pathway missing its last reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


``EnergyMissingHalf`` - Reject an energy metabolism related pathway when half of its non spontaneous reactions lacks a catalyzis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**********Rationale**:

``PathwayScore`` - Accept the pathway if its pathway score exceeds a threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

How is the PathwayScore computed?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :math:`PathwayScore` is a score that aims at estimating how likely a
pathway is to be effectively in the target organism metabolism.

The pathway score is a weighted sum of the :math:`ReactionScore`.

.. math::

   PathwayScore(Pathway) = \frac{\sum_{r \in Reaction(Pathway)} ReactionScore(r)}{|Reaction(Pathway)|}

The (ReactionScore) is computed as follows:

.. math::

    ReactionScore(reaction) = PresenceScore(reaction) +
   KeyReactionScore(reaction) + UniquenessScore(reaction)


The reaction score equals 0 if we did not find any evidence of a
catalysis for this reaction.

Then, the (ReactionScore) components are computed as follows:

.. math::


   PresenceScore(reaction) =
   \begin{cases}
   0.2 & \text{ if we found enzymes catalyzing this reaction } \\
   0 & \text{ otherwise }
   \end{cases}

The key reaction score is computed as

.. math::


   KeyReactionScore(reaction) =
   \begin{cases}
   0.5 & \text{ if the reaction is considered a 'key' reaction for the pathway in MetaCyc } \\
   0 & \text{ otherwise }
   \end{cases}

The uniqueness score is computed as follows:

.. math::


   UniquenessScore(reaction) =  - \exp(|Pathway(reaction)| / 10)

where (\|Pathway(reaction)\|) is the number of pathway having the
reaction (reaction).

How to choose a good threshold?
'''''''''''''''''''''''''''''''

Acknoledgments
--------------

Some of the decision rules implemented in ``pan2met`` are inspired by
the decision rules used by the PathoLogic algorithm.

What could cause an erroneous pathway prediction?
-------------------------------------------------

Unfortunately, there is are a lot of causes that might lead to an
incorrect decision.

An incorrect decision is of two possible kinds: A false positive is a
pathway predicted to be in the metabolism, whereas it is not indeed
realized by the organism metabolism. A false negative, on the other
hand, is a pathway that is rejected by ``pan2met`` decision rules,
whereas it is indeed found in the organism's metabolism.

Incorrect protein catalyzis annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pan2met`` relies on sequence homology method to associate biochemical
reaction catalyzis to a protein. Two protein with similar sequence might
however have distinct biochemical functions, which might fool pan2met
into predicting more pathway than it should.

Gene regulation
~~~~~~~~~~~~~~~

Even if the genome of an organism contains genes coding for enzymes of
all catalyzed steps of a pathway, this does not necessarily implies that
the pathway is indeed effective in the organism. For instance, there
might be gene regulatory mechanism that disable synchronous enzyme
expression thus limiting the potential of such a pathway.
