-- name: insert-pathway(name)$
-- Insert a pathway into the database
replace into pathway (name)
values (:name)
returning id;

-- name: insert-reaction(name)$
-- Insert a reaction into the database
replace into reaction (name)
values (:name)
returning id;

-- name: insert-crossref(local_table, local_id, reference_database, reference_id)^
-- Insert a reference id to table_name(table_id) referencing a record in database identified by reference_soruce, with reference_id
insert into crossref (local_table, local_id, reference_database, reference_id)
values (
    :local_table,
    :local_id,
    :reference_database,
    :reference_id
);

-- name: add-reaction-to-pathway(pathway_id, reaction_id)^
-- Add a reaction to the list of reactions of a given pathway
insert into pathway_reactions (pathway_id, reaction_id)
values (
   :pathway_id,
   :reaction_id
);
