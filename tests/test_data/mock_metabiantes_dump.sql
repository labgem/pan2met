-- Mock (synthetic) data dump compatible with the `metabiantes` schema
-- (see metabiantes/sql/create_schema.sql)
--
-- All identifiers and values are synthetic. No real MetaCyc/BioCyc data is
-- included. Use for tests and local development only.

-- substrate
INSERT INTO substrate (id)
VALUES
 ('CPD-1'),
 ('CPD-2'),
 ('CPD-3'),
 ('CPD-4'),
 ('CPD-5'),
 ('CPD-6');

-- compound
INSERT INTO compound (id, type, comment, atomic_number, atom_charges, smiles, molecular_weight, monoisotopic_mw, gibbs_free_energy)
VALUES
 ('CPD-1', 'Compound', 'Synthetic substrate A', 12, 0, 'CCO', 46.07, 46.0419, -620.0),
 ('CPD-2', 'Compound', 'Synthetic substrate B', 10, -2, 'CC(=O)O', 60.05, 60.0211, -780.0),
 ('CPD-3', 'Compound', 'Synthetic product C', 13, 0, 'CC(=O)OC', 74.08, 74.0368, -590.0),
 ('CPD-4', 'Compound', 'Synthetic carrier', 6, 0, 'CN', 31.06, 31.0422, -120.0),
 ('CPD-5', 'Compound', 'Synthetic energy cofactor', 15, -3, 'OC1C', 150.09, 149.045, -1400.0),
 ('CPD-6', 'Compound', 'Synthetic output product', 8, 0, 'CO', 32.04, 32.0262, -175.0);

-- compound_name_synonymes
INSERT INTO compound_name_synonymes (compound_id, name)
VALUES
 ('CPD-1', 'Substrate A'),
 ('CPD-1', 'Ethanol-like A'),
 ('CPD-2', 'Substrate B'),
 ('CPD-2', 'Acetate-like B'),
 ('CPD-3', 'Product C'),
 ('CPD-4', 'Carrier Molecule'),
 ('CPD-5', 'Energy Cofactor'),
 ('CPD-6', 'Output Product');

-- reaction
INSERT INTO reaction (id, type, comment, ec_number, spontaneous, gibbs_free_energy, physiologically_relevant, reaction_balance_status, reaction_physiological_direction)
VALUES
 ('RXN-0', 'Biochemical-Reaction', 'First step, substrate A -> substrate B', 'EC-1.1.1.1', TRUE, -160.0, TRUE, TRUE, 'left-to-right'),
 ('RXN-1', 'Biochemical-Reaction', 'Second step, substrate B -> product C', 'EC-1.1.1.2', FALSE, 25.0, TRUE, TRUE, 'right-to-left'),
 ('RXN-2', 'Chemical-Reaction', 'Synthetic coupling reaction', NULL, TRUE, -45.0, FALSE, TRUE, 'left-to-right'),
 ('RXN-3', 'Biochemical-Reaction', 'Single substrate reaction', 'EC-2.7.1.1', FALSE, -12.5, TRUE, TRUE, 'left-to-right'),
 ('RXN-4', 'Biochemical-Reaction', 'An example of orphan reaction', 'EC-0.1.2.3', TRUE, -12.5, TRUE, TRUE, 'left-to-right')
 ;

-- reaction_name
INSERT INTO reaction_name (reaction_id, name)
VALUES
 ('RXN-2', 'Synthetic Step 1'),
 ('RXN-1', 'Synthetic Step 2'),
 ('RXN-2', 'Synthetic Coupling'),
 ('RXN-3', 'Synthetic Single-Substrate');

-- reaction_species
INSERT INTO reaction_species (reaction_id, species_id)
VALUES
 ('RXN-0', 1234),
 ('RXN-1', 1234),
 ('RXN-3', 5678);

-- reaction_substrate
INSERT INTO reaction_substrate (reaction_id, substrate_id, stoechiometry, reaction_side)
VALUES
 ('RXN-0', 'CPD-1', 1, 'left'),
 ('RXN-0', 'CPD-5', 1, 'left'),
 ('RXN-0', 'CPD-2', 1, 'right'),
 ('RXN-0', 'CPD-4', 1, 'right'),
 ('RXN-1', 'CPD-2', 1, 'left'),
 ('RXN-1', 'CPD-3', 1, 'right'),
 ('RXN-2', 'CPD-3', 1, 'left'),
 ('RXN-2', 'CPD-6', 1, 'right'),
 ('RXN-3', 'CPD-4', 1, 'left');

-- polypeptide
INSERT INTO polypeptide (id, type, comment, experimental_molecular_weight, molecular_weight, molecular_weight_sequence, half_life, gene, neidhardt_spot_number, atom_charges, isoelectric_point)
VALUES
 ('MONOMER-1', 'monomer', 'Synthetic enzyme subunit alpha', 42000.0, 41800.0, 41950.0, 600.0, 'synA', NULL, -10, 5.1),
 ('MONOMER-2', 'monomer', 'Synthetic enzyme subunit beta', 35000.0, 34800.0, 34900.0, 900.0, 'synB', NULL, -6, 6.4),
 ('CPLX-1', 'complex', 'Synthetic protein complex', 77000.0, 76600.0, NULL, NULL, NULL, NULL, NULL, NULL);

-- polypeptide_complex_component
INSERT INTO polypeptide_complex_component (complex_id, component_id, coefficient)
VALUES
 ('CPLX-1', 'MONOMER-1', 1),
 ('CPLX-1', 'MONOMER-3', 1);

-- reaction_enzyme
INSERT INTO reaction_enzyme (reaction_id, enzyme_id)
VALUES
 ('RXN-0', 'MONOMER-1'),
 ('RXN-1', 'MONOMER-2'),
 ('RXN-2', 'CPLX-1');

-- pathway
INSERT INTO pathway (id, comment)
VALUES
 ('PWY-0', 'Synthetic core metabolic pathway'),
 ('PWY-1', 'Synthetic downstream pathway'),
 ('PWY-2', 'Synthetic variant pathway');

-- pathway_name_synonymes
INSERT INTO pathway_name_synonymes (pathway_id, name)
VALUES
 ('PWY-0', 'Synthetic Core'),
 ('PWY-1', 'Synthetic Downstream'),
 ('PWY-2', 'Synthetic Variant');

-- pathway_sub_pathway
INSERT INTO pathway_sub_pathway (super_pathway_id, sub_pathway_id)
VALUES
 ('PWY-1', 'PWY-2');

-- pathway_reaction_graph
INSERT INTO pathway_reaction_graph (pathway_id, predecessor_reaction_id, successor_reaction_id)
VALUES
 ('PWY-0', 'RXN-0', 'RXN-1'),
 ('PWY-0', 'RXN-1', 'RXN-2');

-- pathway_reaction
INSERT INTO pathway_reaction (pathway_id, reaction_id, reaction_direction)
VALUES
 ('PWY-0', 'RXN-0', 'left-to-right'),
 ('PWY-0', 'RXN-1', 'right-to-left'),
 ('PWY-0', 'RXN-2', 'left-to-right'),
 ('PWY-1', 'RXN-3', 'left-to-right'),
 ('PWY-1', 'RXN-0', 'left-to-right')
 ;

-- pathway_key_reaction
INSERT INTO pathway_key_reaction (pathway_id, reaction_id)
VALUES
 ('PWY-0', 'RXN-0'),
 ('PWY-1', 'RXN-3');

-- pathway_species
INSERT INTO pathway_species (pathway_id, species_id)
VALUES
 ('PWY-0', 1234),
 ('PWY-1', 12345);

-- pathway_taxonomic_range
INSERT INTO pathway_taxonomic_range (pathway_id, taxon_id)
VALUES
 ('PWY-0', 1234),
 ('PWY-1', 5678);

-- pathway_variant
INSERT INTO pathway_variant (pathway_id, variant_id)
VALUES
 ('PWY-0', 'PWY-2');

-- pathway_variant_group
INSERT INTO pathway_variant_group (variant_group_id, variant_id)
VALUES
 ('1', 'PWY-0'),
 ('1', 'PWY-2');

-- pathway_ontology
INSERT INTO pathway_ontology (pathway_id, path, depth, pathway_class)
VALUES
 ('PWY-0', 1, 0, 'Root'),
 ('PWY-0', 2, 1, 'Synthetic'),
 ('PWY-0', 3, 2, 'SyntheticMetabolism'),
 ('PWY-1', 1, 0, 'Root'),
 ('PWY-1', 2, 1, 'Synthetic');
