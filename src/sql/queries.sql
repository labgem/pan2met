-- name: get_all_pathways()
-- Get all metabolic pathways from the database
SELECT DISTINCT subject
FROM statements
WHERE object = 'bp:Pathway';

-- name: get_biochemical_reactions_by_pathway
-- Get all reactions of a metabolic pathway
SELECT DISTINCT object
FROM statements
WHERE predicate = 'bp:pathwayComponent'
AND subject = :pathway
AND object LIKE 'bp:BiochemicalReaction%';

-- name: get_biochemical_reactions_id_by_pathway
-- Get all BioPAX pb:id of a reaction of a pathway
SELECT label.object, value
FROM statements label
WHERE predicate = 'bp:id'
AND label.stanza IN
(
SELECT DISTINCT reaction.object
FROM statements reaction
WHERE reaction.predicate = 'bp:pathwayComponent'
AND reaction.subject = :pathway
AND reaction.object LIKE 'bp:BiochemicalReaction%');

-- name: get_reaction_id_by_reaction(reaction)
-- Get the first ID associated with the reaction object (e.g., bp:BiochemicalReactionXXXX -> RXN-1.2.3.1)
SELECT id.value, reaction.subject
FROM statements id
INNER JOIN statements xref
ON xref.subject = id.subject
INNER JOIN statements reaction
ON reaction.subject = xref.subject
WHERE id.predicate = 'bp:id'
AND xref.predicate = 'bp:xref'
AND reaction.subject = :reaction;


SELECT id.value
FROM statements id
INNER JOIN statements id_type
ON id_type.subject = id.subject
WHERE id_type.predicate = 'rdf:type'
AND id_type.object = 'bp:UnificationXref';

SELECT id.value
FROM statements id
INNER JOIN statements id_type
ON id_type.subject = id.subject
INNER JOIN statements reaction
ON
WHERE id_type.predicate = 'rdf:type'
AND id_type.object = 'bp:UnificationXref';
