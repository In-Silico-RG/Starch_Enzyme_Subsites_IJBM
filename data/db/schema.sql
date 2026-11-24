CREATE TABLE bond (pdb, conn_type, instance_id_1, chain_1, seq_1, comp_1, atom_1, instance_id_2, chain_2, seq_2, comp_2, atom_2, distance_A);

CREATE TABLE census_linkage ('pdb', 'family', 'EC', 'protein', 'organism', 'method', 'resolution_A', 'nonreducing_ring_C1', 'reducing_ring_O4', 'phi_deg', 'psi_deg', 'theta_nonred_deg', 'theta_red_deg', 'conformation', 'rscc_nonred', 'rscc_red', 'occ_nonred', 'occ_red', 'density', 'protein_atoms_within_4A', 'O4_to_nearest_carboxylate_A', 'nearest_carboxylate');

CREATE TABLE census_on_maps ('pdb', 'family', 'nonreducing_ring_C1', 'reducing_ring_O4', 'phi_deg', 'psi_deg', 'conformation', 'density', 'theta_nonred_deg', 'theta_red_deg', 'resolution_A', 'G_CHARMM36_kJ', 'G_GLYCAM06_kJ', 'E_xtb_kJ');

CREATE TABLE chembl_activity (activity_id, target_chembl_id, molecule_chembl_id, canonical_smiles,
                   standard_type, standard_relation, standard_value, standard_units, pchembl_value, assay_chembl_id,
                   assay_type, assay_description, document_chembl_id, document_year, data_validity_comment, target_organism);

CREATE TABLE chembl_molecule (molecule_chembl_id PRIMARY KEY, pref_name, canonical_smiles, inchikey,
                   full_mwt, alogp, hbd, hba, num_ro5_violations, max_phase);

CREATE TABLE chembl_pdb_ligand(
  comp_id,
  pdb_name,
  molecule_chembl_id,
  pref_name
);

CREATE TABLE chembl_target (uniprot, target_chembl_id, target_type, pref_name, organism);

CREATE TABLE component (comp_id PRIMARY KEY, name, type, formula, formula_weight, smiles, inchikey);

CREATE TABLE contact (instance_id, partner_kind, partner_instance_id, chain, seq, ins_code, comp, min_distance_A, atom, partner_atom);

CREATE TABLE contact_uniprot(
  contact_rowid INT,
  instance_id,
  pdb,
  chain,
  seq,
  ins_code,
  comp,
  uniprot,
  uniprot_resnum
);

CREATE TABLE entry (pdb PRIMARY KEY, method, resolution_A, release_date, family, EC, protein, organism, uniprot, n_waters);

CREATE TABLE instance (instance_id PRIMARY KEY, pdb, auth_asym, label_asym, auth_seq, ins_code, comp_id, entity_id, kind, n_atoms, altlocs, mean_occupancy, mean_b, rscc, rsr, x, y, z, ligand_class);

CREATE TABLE ligand_entity ('pdb', 'family', 'EC', 'protein', 'organism', 'kind', 'ligand_class', 'entity', 'comp_ids', 'n_residues', 'name', 'prd', 'links', 'formula', 'smiles', 'wurcs', 'copies', 'chains');

CREATE TABLE part1_activity ('group', 'family', 'EC', 'name', 'attack', 'bond', 'reactant', 'reason', 'mechanism', 'clan', 'fold', 'nucleophile_base', 'proton_donor', 'n_sequences');

CREATE TABLE part1_excluded_activity ('family', 'EC', 'name', 'bond', 'reason');

CREATE TABLE part2_pdb_detail ('pdb', 'method', 'release_date', 'resolution_A', 'longest_glucan_residues', 'branched_oligosaccharides', 'other_ligands');

CREATE TABLE part2_structure ('family', 'in_scope', 'EC', 'protein', 'organism', 'uniprot', 'pdb', 'chains', 'carbohydrate_ligands', 'resolution');

CREATE TABLE part3_catalytic_residue ('family', 'mechanism_CAZy', 'nucleophile_or_base_CAZy', 'acid_CAZy', 'reference_enzyme', 'uniprot', 'catalytic_residues_reference_numbering', 'residue_source', 'candidate_reference_structures_PDB,res,longest_bound_glucan', 'notes');

CREATE TABLE part3_site_geometry ('family', 'structure', 'pdb', 'mechanism', 'residue_source', 'nuc_or_base', 'acid', 'd_CC_A', 'd_OO_min_A', 'ligand_at__1', 'd_nuc_base_C_to_C1_A', 'plus1', 'd_acid_O_to_X_A', 'Q_A', 'theta_deg', 'phiCP_deg', 'pucker_at__1', 'phi_deg', 'psi_deg', 'linkage_note');

CREATE TABLE sifts_domain (pdb, chain, db, accession, identifier, auth_start, auth_end, domain);

CREATE TABLE sifts_residue (pdb, chain, auth_seq, ins_code, resname, uniprot, uniprot_resnum, uniprot_res);

CREATE INDEX i_bond_pdb ON bond(pdb);

CREATE INDEX i_cont_inst ON contact(instance_id);

CREATE INDEX i_inst_comp ON instance(comp_id);

CREATE INDEX i_inst_pdb ON instance(pdb);

CREATE INDEX i_sres ON sifts_residue(pdb, chain, auth_seq, ins_code);

