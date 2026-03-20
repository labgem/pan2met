-- name: get_pathways()
-- Get all pathway in the database
SELECT name FROM pathway;

-- name: get_reactions()
-- Get all reactions in the database
SELECT name FROM reaction;

-- name: get_orphan_reactions()
-- Get all orphan reactions in the database
SELECT name FROM reaction
WHERE NOT EXISTS (
    SELECT 1
    FROM reaction_enzyme
    WHERE reaction_enzyme.reaction_id = reaction.id
);


-- name: get_spontaneous_reactions()
-- Get all spontaneous reactions in the database
SELECT name FROM reaction
WHERE spontaneous = TRUE;


-- name: get_reactions_of_pathway(pathway_id)
-- Get all reactions of a pathway
SELECT reaction.name
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE pathway.name = :pathway_id;

-- name: get_key_reactions_of_pathway(pathway_id)
-- Get all reactions of a pathway
SELECT reaction.name
FROM (SELECT pathway.id, pathway.name
        FROM pathway
        WHERE name = :pathway_id)
        pathway
INNER JOIN pathway_key_reaction
ON pathway_key_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_key_reaction.reaction_id;

-- name: get_non_spontaneous_reactions_of_pathway(pathway_id)
-- Get all non-spontaneous reactions of a pathway
SELECT reaction.name
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE pathway.name = :pathway_id
WHERE reaction.spontaneous IS DISTINCT FROM TRUE;

-- name: get_non_orphan_non_spontaneous_reactions_of_pathway(pathway_id)
-- Get all non-spontaneous and non-orphan reactions of a pathway
SELECT reaction.name
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE pathway.name = :pathway_id
AND (reaction.spontaneous IS DISTINCT FROM TRUE
OR NOT EXISTS (
    SELECT 1
    FROM reaction_enzyme
    WHERE reaction_enzyme.reaction_id = reaction.id
));

-- name: get_enzymes_of_reaction(reaction_id)
-- List the enzymes catalyzing a reaction
SELECT polypeptide.name
FROM reaction
INNER JOIN reaction_enzyme
ON reaction_enzyme.reaction_id = reaction.reaction_id
INNER JOIN polypeptide
ON polypeptide.id = reaction_enzyme.enzyme_id
WHERE reaction.name = :reaction_id;

-- name: get_reactions_by_ec_number(ec_number)
-- List the reactions having the given EC number
SELECT reaction.name
FROM reaction
WHERE reaction.ec_number = :ec_number;

-- name: get_pathway_taxonomic_range(pathway_id)^
-- Get the pathway taxonomic range NCBI-Taxonomy identifier
SELECT pathway_taxonomic_range.taxon_id
FROM pathway
INNER JOIN pathway_taxonomic_range
ON pathway_taxonomic_range.pathway_id = pathway.id
WHERE pathway.name = :pathway_id;

-- name: reaction_is_key(pathway_id, reaction_id)^
-- Check if the reaction is a key for the pathway
SELECT 1
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE pathway.name = :pathway_id
AND reaction.name = :reaction_id;

-- name: get_pathways_with_reaction(reaction_id)
-- List all pathways
SELECT pathway.name
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE reaction.name = :reaction_id;

-- name: count_pathways_with_reaction(reaction_id)^
-- List all pathways
SELECT count(pathway.name)
FROM pathway
INNER JOIN pathway_reaction
ON pathway_reaction.pathway_id = pathway.id
INNER JOIN reaction
ON reaction.id = pathway_reaction.reaction_id
WHERE reaction.name = :reaction_id;

-- name: get_variants_of_pathway(pathway_id)
-- List all variants of a pathway
SELECT variant.name
FROM pathway, pathway variant, pathway_variant
WHERE pathway_variant.pathway_id = pathway.id
AND pathway_variant.variant_id = variant.id
AND pathway.name = :pathway_id;

-- name: get_ontology_parent_class_of_pathway(pathway_id)
-- List all ontology parent class of a pathway
SELECT pathway_ontology.pathway_class
FROM pathway_ontology
INNER JOIN pathway
ON pathway.id = pathway_ontology.pathway_id
WHERE pathway.name = :pathway_id;


-- name: get_pathway_reaction_order(pathway_id)
-- List all pairs of reactions in the pathway, where the first reaction is a predecessor of the
-- second reaction in the pathway.
SELECT predecessor_reaction.name, successor_reaction.name
FROM pathway
INNER JOIN pathway_reaction_graph
ON pathway_reaction_graph.pathway_id = pathway.id
INNER JOIN reaction AS predecessor_reaction
ON predecessor_reaction.id = pathway_reaction_graph.predecessor_reaction_id
INNER JOIN reaction AS successor_reaction
ON successor_reaction.id = pathway_reaction_graph.successor_reaction_id
WHERE pathway.name = :pathway_id;
