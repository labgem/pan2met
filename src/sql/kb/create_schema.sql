-- name: create_schema#
-- Prepare the MetaCyc SQL database schema
create table protein_reference (
    id integer not null,
    sequence text not null,
    constraint pk_protein_reference primary key (id)
);

create table protein_monomer (
    id integer not null,
    protein_reference_id integer,
    constraint pk_protein_monomer primary key (id),
);

create table protein_complex_member (
   complex_id integer not null,
   protein_monomer_id integer not null,
   constraint fk_protein_complex_member foreign key protein_monomer_id references protein_monomer (id),
   constraint primary key (complex_id, protein_monomer_id)
);

create table reaction (
    id integer not null,
    name text,
    ec_number text,
    is_spontaneous integer, -- bool (0, 1)
    is_orphan integer, -- bool (0, 1)
    constraint pk_reaction primary key (id),
    constraint unique_name_reaction unique (name)
);

create table pathway (
    id integer not null,
    name text,
    constraint pk_pathway primary key (id),
    constraint unique_name_pathway unique (name)
);

create table pathway_reactions (
    pathway_id not null,
    reaction_id not null,
    foreign key (pathway_id) references pathway(id),
    foreign key (reaction_id) references reaction(id),
    constraint pk_pathway_reactions primary key (pathway_id, reaction_id)
);

create table crossref (
    id integer not null,
    local_table text not null,
    local_id integer not null,
    reference_database text not null,
    reference_id text not null,
    constraint pk_crossref primary key (id)
);
