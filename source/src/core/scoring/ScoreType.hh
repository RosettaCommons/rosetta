// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreType.hh
/// @brief  Score type enumeration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_ScoreType_hh
#define INCLUDED_core_scoring_ScoreType_hh

#include <iostream>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

	/////////////////////////////////////////////////////////////////////////////////
	/////// WARNING WARNING WARNING
	///////
	/////// if you add a new ScoreType please also add its string name in ScoreTypeManager.cc
	///////
	/////// WARNING WARNING WARNING
	/////////////////////////////////////////////////////////////////////////////////

/// @brief Type for looking up cached energies
/// I guess we could get rid of the fa_ prefix, except maybe for
/// fa_pair, to distinguish from std::pair and the centroid pair score...
enum ScoreType {

	/// begin short ranged ci2b scores -- these guys are cached
	/// in the energy graph -- when appropriate --
	/// they are reused between rounds of scoring.
	fa_atr = 1, //enumeration starts at 1 for indexing utility::vector1
	fa_rep,
	fa_sol,
	fa_intra_atr,
	fa_intra_rep,
	fa_intra_sol,
	fa_intra_RNA_base_phos_atr, //RNA specific score term
	fa_intra_RNA_base_phos_rep, //RNA specific score term
	fa_intra_RNA_base_phos_sol, //RNA specific score term
	lk_hack,
	lk_ball,
	lk_ball_wtd,
	lk_ball_iso,
	coarse_fa_atr,
	coarse_fa_rep,
	coarse_fa_sol,
	coarse_beadlj,
	mm_lj_intra_rep,
	mm_lj_intra_atr,
	mm_lj_inter_rep,
	mm_lj_inter_atr,
	mm_twist,      // could be lr 2benergy and not in energy graph
	mm_bend,       // could be lr 2benergy and not in energy graph
	mm_stretch,    // could be lr 2benergy and not in energy graph
	lk_costheta,
	lk_polar,
	lk_nonpolar,
	lk_polar_intra_RNA,    //RNA specific score term
	lk_nonpolar_intra_RNA, //RNA specific score term
//	csd_torsion, //commenting out until it is implemented
	fa_elec,
	fa_elec_bb_bb,
	fa_elec_bb_sc,
	fa_elec_sc_sc,
	h2o_hbond,
	dna_dr,
	dna_bp,
	dna_bs,
	peptide_bond,
	pcs, //Pseudocontact Shift Energy
	pcs2, //Pseudocontact Shift Energy version 2. Will replace pcs end of 2010

	fastsaxs,   // fastsaxs agreement using formulation of Stovgaard et al (BMC Bioinf. 2010)
	saxs_score, // centroid saxs asessment
	saxs_cen_score,
	saxs_fa_score, // full-atom SAXS score
	pddf_score, // score based on pairwise distance distribution function

	//pba Membrane all atom terms
	fa_mbenv,       // depth dependent reference term
	fa_mbsolv,      // burial+depth dependent term

	//Split out fa_elec for RNA.
	fa_elec_rna_phos_phos,
	fa_elec_rna_phos_sugr,
	fa_elec_rna_phos_base,
	fa_elec_rna_sugr_sugr,
	fa_elec_rna_sugr_base,
	fa_elec_rna_base_base,

	fa_elec_aro_aro,
	fa_elec_aro_all,
	hack_aro,

	rna_fa_atr_base,
	rna_fa_rep_base,
	rna_data_backbone, // Using chemical accessibility data for RNA.

	//Carbon hydrogen bonds -- moved up here since they are currently context independent.
	ch_bond,
	ch_bond_bb_bb,
	ch_bond_sc_sc,
	ch_bond_bb_sc,

	/// proline closure energy
	pro_close,
	rama2b,

	vdw,     // centroid
	cenpack, // centroid

	cenpack_smooth,   //fpd  smooth cenpack
	cen_hb,           //fpd  centroid bb hbonding

	hybrid_vdw,     // hybrid centroid+fa

	// gaussian overlap
	gauss,


	rna_vdw,          // low res clash check for RNA
	rna_base_backbone,          // Bases to 2'-OH, phosphates, etc.
	rna_backbone_backbone,      // 2'-OH to 2'-OH, phosphates, etc.
	rna_repulsive,              // mainly phosphate-phosphate repulsion
	rna_base_pair_pairwise,    // Base-base interactions (Watson-Crick and non-Watson-Crick)
	rna_base_axis_pairwise,    // Force base normals to be parallel
	rna_base_stagger_pairwise, // Force base pairs to be in same plane.
	rna_base_stack_pairwise,   // Stacking interactions
	rna_base_stack_axis_pairwise,   // Stacking interactions should involve parallel bases.
	rna_data_base, // Using chemical accessibility data for RNA.

	//RNA stuff
	// This is a filtered version of the pairwise RNA low resolution terms above,
	//  disallows a base edge to form  more than one base pair, and
	//  disallows two bases to both stack and pair.
	// THIS IS NOT REALLY PAIR-WISE, but is calculated in a finalize_energy
	//  step at the end of a 2-body  score function (RNA_PairwiseLowResolutionEnergy.cc)
	rna_base_pair,    // Base-base interactions (Watson-Crick and non-Watson-Crick)
	rna_base_axis,    // Force base normals to be parallel
	rna_base_stagger, // Force base pairs to be in same plane.
	rna_base_stack,   // Stacking interactions
	rna_base_stack_axis,   // Stacking interactions should involve parallel bases.

  rna_mg, //knowledge-based term for mg(2+)/RNA interactions for use in low res modeling.
  rna_mg_rep, // ad-hoc, empirically validated term to prevent uncommon mg(2+)/atom interactions.
  rna_mg_indirect, //knowledge-based term for mg(2+)/RNA interactions for use in low res modeling.

	// High resolution
	rna_torsion,       // RNA torsional potential.
	rna_sugar_close,   // constraints to keep RNA sugar closed, and with reasonably ideal geometry
	fa_stack,          // stacking interaction modeled as pairwise atom-atom interactions
	fa_stack_aro,
	stack_elec,          // distance dependent dielectric between base atoms (attenuated parallel to plane)
	stack_elec_base_base,
	stack_elec_base_bb,

	// DNA constraints-based torsional potentials
	dna_bb_torsion,
	dna_sugar_close,
	dna_base_distance,

	geom_sol_fast,           //Context independent version. Currently tested only for RNA case.
	geom_sol_fast_intra_RNA, //RNA specific score term

	fa_cust_pair_dist,  // custom short range 2b
	custom_atom_pair,


	// All the orbitals scoretypes
	orbitals_hpol_bb,
	pci_cation_pi,
	pci_pi_pi,
	pci_salt_bridge,
	pci_hbond,


	#ifdef PYROSETTA
		PyRosettaTwoBodyContextIndepenedentEnergy_first,
		PyRosettaTwoBodyContextIndepenedentEnergy_last = PyRosettaTwoBodyContextIndepenedentEnergy_first + 10,
	#endif

	// for in-python runtime-defined methods: most inclusive positioning in the
	// score type enum: user may define any kind of energy method
	// However, the user may define only a single runtime defined method
	python,  // <-- Deprecated use PyRosettaEnergie* instead

	n_ci_2b_score_types = python, /// keep this guy at the end of the ci2b scores
	//end ci2b scores

	// Begin short-ranged, context dependent two-body energy method types.
	// These are also cached in the edges of the EnergyGraph.
	fa_pair, /// == fa_pair_pol_pol
	fa_pair_aro_aro,
	fa_pair_aro_pol,
	fa_pair_pol_pol,
	fa_plane,
	hbond_sr_bb,
	hbond_lr_bb,
	hbond_bb_sc,
	hbond_sr_bb_sc,
	hbond_lr_bb_sc,
	hbond_sc,
	hbond_intra,           //Currently effects only RNA

	#ifdef PYROSETTA
		PyRosettaTwoBodyContextDependentEnergy_first,
		PyRosettaTwoBodyContextDependentEnergy_last = PyRosettaTwoBodyContextDependentEnergy_first + 10,
	#endif

	// protein-protein interface scores
	interface_dd_pair,

	// Geometric solvation
	geom_sol,
	geom_sol_intra_RNA,    //RNA specific score term
	occ_sol_fitted,
	occ_sol_fitted_onebody,
	occ_sol_exact,

	// centroid rotamer pair, P(r,ang,dih|aa)
	cen_rot_pair, //P(r|aa)
	cen_rot_pair_ang, //P(ang|r,aa)
	cen_rot_pair_dih, //P(dih|r,aa)
	pair, // centroid
	cen_pair_smooth,  //fpd  smooth centroid pair
	Mpair,
	// sucker atom energy
	suck,

	//RNA stuff
	//Low resolution
	rna_rg,           // Radius of gyration for RNA

	// nucleotide resolution thermodynamics
	rna_loop,  // RNA loop closure terms -- attempting model full RNA folding free energy
	rna_loop_fixed,
	rna_loop_logN,
	rna_loop_harmonic,

	//  FACTS solvation model
	facts_elec,
	facts_solv,
	facts_sasa,

	// centroid interchain 1b (docking) scores
	interchain_pair,
	interchain_vdw,

	//
	// end short ranged two body scores
	n_shortranged_2b_score_types = interchain_vdw, // keep this guy at the end of the sr ci/cd 2b scores
	// 30 as of 1/7/2007 -- don't ever code using the literal "30", this is just a helpful count

	gb_elec,

	//Full atom disulfide terms
	dslf_ss_dst,
	dslf_cs_ang,
	dslf_ss_dih,
	dslf_ca_dih,
	dslf_cbs_ds,
	// supercedes the above
	dslf_fa13,
	//Centroid disulfide terms
	dslfc_cen_dst,
	dslfc_cb_dst,
	dslfc_ang,
	dslfc_cb_dih,
	dslfc_bb_dih,
	//disulfide matching terms
	dslfc_rot,
	dslfc_trans,
	dslfc_RT,

	atom_pair_constraint,
	constant_constraint,
	coordinate_constraint,
	angle_constraint,
	dihedral_constraint,
	big_bin_constraint,
	dunbrack_constraint,
	site_constraint,
	metalhash_constraint, // Rigid body, metal binding constraints for centroid mode
	rna_bond_geometry, // deviations from ideal geometry

	rama,
	omega,
	fa_dun,
	fa_dun_dev,
	fa_dun_rot,
	fa_dun_semi,
	dna_chi,
	p_aa_pp,
	p_aa_pp_offset,
	yhh_planarity,
	h2o_intra,
	ref,
	ref_nc,
	seqdep_ref,
  nmer_ref,
  nmer_pssm,
  nmer_svm,
	envsmooth,
	e_pH,
	rna_bulge,
  mg_ref,  // chemical potential for mg(2+) ('reference weight' in Rosetta lingo)
	free_P, // bonus for virtualizing RNA phosphate
	free_2HOstar, // bonus for virtualizing RNA 2'-OH
	intermol, // cost of instantiating a chain form 1 M std state.
	special_rot,

	other_subpose, // in preparation for multi-pose stuff.

	// PB potential
	PB_elec,

	/// Whole structure energies
	/// centroid whole structure energies
	cen_env_smooth,   //fpd smooth centroid env
	cbeta_smooth,     //fpd smooth cbeta
	cen_rot_env,
	cen_rot_dun,
	env,
	cbeta,
	DFIRE,
	Menv,
	Mcbeta,
	Menv_non_helix,
	Menv_termini,
	Menv_tm_proj,
	Mlipo,
	rg, // radius of gyration
	rg_local, //radius of gyration for repeat proteins
	co, // contact order
	hs_pair,
	ss_pair,
	rsigma,
	sheet,
	burial, // informatic burial prediction
	abego,  // informatic torsion-bin prediction

	/// Whole structure energies, centroid score
	// secondary structure scores
	natbias_ss,
	natbias_hs,
	natbias_hh,
	natbias_stwist,

	/// amino acid composition score
	aa_cmp,

	dock_ens_conf, //conformer reference energies for docking

	csa,//NMR chemical shift anisotropy energy
	dc,//NMR dipolar coupling energy
	rdc,//NMR residual dipolar coupling energy
	rdc_segments, //fit alignment on multiple segments independently
	rdc_rohl,
  // end centroid whole structure energies
	holes,
	holes_decoy,
	holes_resl,
	holes_min,
	holes_min_mean,

	rna_chem_shift, //RNA NMR chemical shift pseudo-energy term

	dab_sasa, // classic 1.4A probe solvant accessible surface area
	dab_sev,  // solvent excluded volume -- volume of atoms inflated by 1.4A
	sa, // nonpolar contribution in GBSA
	d2h_sa, // correlation between sasa and hydrogen exchange data.
	ProQM, //membrane MQAP.
	ProQ, // MQAP.

	// centroid interhcain 1b (docking) scores -- Demonstrate that these improve performance or remove them
	interchain_env,
	interchain_contact,
	//

	chainbreak,
	linear_chainbreak,
	overlap_chainbreak,
	distance_chainbreak,
	dof_constraint,

	rama2b_offset,
	omega2b_offset,

	cart_bonded,  // cartesian bonded potential
	cart_bonded_angle,  // cartesian bonded potential
	cart_bonded_length,  // cartesian bonded potential
	cart_bonded_torsion,  // cartesian bonded potential

	//Neighbor Vector solvation approximation
	neigh_vect,
	neigh_count,
	neigh_vect_raw,

	//Symmetry bonus
	symE_bonus,

	//Implicit Ligand interactions (symmetry)
	sym_lig,

	// Other energies.

	// packing score energy
	pack_stat,

	// model-quality metrics.
	rms,


	// for ResidueConstraint
	res_type_constraint,

  // Residue Type linking constraint
	res_type_linking_constraint,

	// for PocketConstraint
	pocket_constraint,

	// for BackboneStubConstraint
	backbone_stub_constraint,
  backbone_stub_linear_constraint,

	surface,
	p_aa,
	unfolded,

	// fit-to-density scores
	elec_dens_fast,
	elec_dens_window,
	elec_dens_whole_structure_ca,
	elec_dens_whole_structure_allatom,
	elec_dens_atomwise,

	// patterson correlation
	patterson_cc,

	// crystallographic ML target
	xtal_ml,

	hpatch,
  //membrane environment smooth
	Menv_smooth,

  #ifdef PYROSETTA
		PyRosettaEnergy_first,
		PyRosettaEnergy_last = PyRosettaEnergy_first + 10,
	#endif

	// etc etc
	// Why is there a total score?
	total_score,

	// Dummy score type to insure that PyRosetta can correctly identify when total_score is used, see bug [bug #0000091] for details
	dummy_score_type,

	/// This element marks the end of the active score types.  Elements in the enumeration
	/// up to this point will have space allocated for them in the EnergyMap object.  Elements
	/// past this point are considered inactive and will not have space allocated for them.
	/// If you wish to use an inactive score type, you must move that score type into its appropriate
	/// position in the ScoreType enumeration (described above) and then recompile.
	/// Inactive score types must still have their names included in the ScoreTypeManager's
	/// string-to-score-type map.
	n_score_types = dummy_score_type,

	/// This element marks the very end of the score type enumeration.  Elements between
	/// the n_score_types element and this element are considered inactive.  They may not
	/// be used by any EnergyMethod or they will result in an out-of-bounds write and
	/// unexpected behavior.  To use an inactived score type, the score type must be moved
	/// to an earlier position in this enumeration, and the program must be recompiled.
	/// Keep this guy last.
	end_of_score_type_enumeration = dummy_score_type
};


typedef utility::vector1< ScoreType > ScoreTypes;

/// @brief Returns the ScoreType titled  <name>
///
/// example(s):
///     score_type_from_name("fa_sol")
/// See also:
///     ScoreFunction
///     ScoreType
///     Energies
///     Energies.residue_total_energies
///     name_from_score_type
ScoreType
score_type_from_name( std::string const & name );

/// @brief Returns the name of the ScoreType  <score_type>
///
/// example(s):
///     name_from_score_type(fa_sol)
/// See also:
///     ScoreFunction
///     ScoreType
///     Energies
///     Energies.residue_total_energies
///     score_type_from_name
std::string
name_from_score_type( ScoreType  score_type );

	/// @brief input operator for ScoreType enum type
std::istream & operator >>( std::istream & is, ScoreType & t );
	/// @brief output operator for ScoreType enum type
std::ostream & operator <<( std::ostream & os, ScoreType const & t );

	/// @brief output operator for ScoreTypes list
std::ostream & operator <<( std::ostream & os, ScoreTypes const & score_types );

} // namespace scoring
} // namespace core


#endif
