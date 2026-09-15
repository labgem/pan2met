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
You can use the file
```conf/example.yaml`` <https://github.com/labgem/pan2met-rs/blob/main/conf/example.yaml>`__

For instance, a configuration file ``configuration.yaml`` with the
following content would enable only the filter on transport pathwways
and signaling pathway (``TransportPathway`` and ``SignalingPathway``
rule identifiers), the filter of pathway with no catalysis evidence
found (``AllReactionsMissing``) and would keep all pathways with all
reactions catalyzed, or spontaneous (``AllReactionsCatalyzed``).

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

``TransportPathway`` - Reject the transport pathways
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SignalingPathway`` - Reject the signaling pathways
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AllReactionsCatalyzed`` - Accept a pathway when all its non spontaneous reactions are catalyzed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AllReactionsMissing`` - Reject a pathway when none of its non spontaneous reactions is catalyzed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``KeyReaction``
^^^^^^^^^^^^^^^

``SynthesisMissingLast`` - Reject a biosynthesis pathway missing its last reaction.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: If the last reaction leading to the metabolite of
interest is missing, it might be a sign that the metabolite is not
synthesized in the organism, hence that the biosynthesis pathway is not
there.

``DegradationMissingFirst`` - Reject a degradation pathway missing its last reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**: If the first reaction in a degradation pathway, the one
that initially transformed the metabolite targetted by the degradation
pathway, is not catalyzed, it might be a sign that no such degradation
occurs in the organism.

``EnergyMissingHalf`` - Reject an energy metabolism related pathway when half of its non spontaneous reactions lacks a catalyzis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Rationale**:

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
