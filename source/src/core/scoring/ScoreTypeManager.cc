// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreTypeManager.cc
/// @brief  Method definitions for ScoreTypeManager
/// @author Andrew Leaver-Fay

// Unit Heaaders
#include <core/scoring/ScoreTypeManager.hh>


// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <map>
#include <string>
#include <iostream>

#include <sstream>
#include <utility/vector1_bool.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {

bool ScoreTypeManager::initialized_( false );
std::map< std::string, ScoreType > ScoreTypeManager::name2score_type_;
utility::vector1< std::string > ScoreTypeManager::score_type2name_;


void fill_score_range(std::map< std::string, ScoreType > & M, std::string prefix, int first, int last)
{
	M[ prefix + "_first" ] = ScoreType(first);
	M[ prefix + "_last" ] = ScoreType(last);
	for ( int i=first+1; i<last; i++ ) {
		std::ostringstream s; s << prefix << '_' << i;
		M[ s.str() ] = ScoreType(i);
	}
}


/// @brief initialize the ScoreType name vector and map
///
/// @details initialize all the SCORETYPE string name into the vector then set up
/// the look-up map from string name to enum type
void
ScoreTypeManager::setup_score_type_names()
{
	if ( initialized_ ) return;
	initialized_ = true;

	name2score_type_[ "fa_atr" ] = fa_atr;
	name2score_type_[ "fa_rep" ] = fa_rep;
	name2score_type_[ "fa_sol" ] = fa_sol;
	name2score_type_[ "lk_hack" ] = lk_hack;
	name2score_type_[ "lk_ball" ] = lk_ball;
	name2score_type_[ "lk_ball_wtd" ] = lk_ball_wtd;
	name2score_type_[ "lk_ball_iso" ] = lk_ball_iso;
	name2score_type_[ "lk_costheta" ] = lk_costheta;
	name2score_type_[ "lk_polar" ] = lk_polar;
	name2score_type_[ "lk_nonpolar" ] = lk_nonpolar;
	name2score_type_[ "fa_intra_atr" ] = fa_intra_atr;
	name2score_type_[ "fa_intra_rep" ] = fa_intra_rep;
	name2score_type_[ "fa_intra_sol" ] = fa_intra_sol;
	name2score_type_[ "fa_intra_atr_xover4" ] = fa_intra_atr_xover4;
	name2score_type_[ "fa_intra_rep_xover4" ] = fa_intra_rep_xover4;
	name2score_type_[ "fa_intra_sol_xover4" ] = fa_intra_sol_xover4;
	name2score_type_[ "fa_atr_dummy" ] = fa_atr_dummy;
	name2score_type_[ "fa_rep_dummy" ] = fa_rep_dummy;
	name2score_type_[ "fa_sol_dummy" ] = fa_sol_dummy;
	name2score_type_[ "fa_vdw_tinker" ] = fa_vdw_tinker;
	name2score_type_[ "coarse_fa_atr" ] = coarse_fa_atr;
	name2score_type_[ "coarse_fa_rep" ] = coarse_fa_rep;
	name2score_type_[ "coarse_fa_sol" ] = coarse_fa_sol;
	name2score_type_[ "coarse_beadlj" ] = coarse_beadlj;
	name2score_type_[ "mm_lj_intra_rep" ] = mm_lj_intra_rep;
	name2score_type_[ "mm_lj_intra_atr" ] = mm_lj_intra_atr;
	name2score_type_[ "mm_lj_inter_rep" ] = mm_lj_inter_rep;
	name2score_type_[ "mm_lj_inter_atr" ] = mm_lj_inter_atr;
	name2score_type_[ "mm_twist" ] = mm_twist;
	name2score_type_[ "mm_bend" ] = mm_bend;
	name2score_type_[ "mm_stretch" ] = mm_stretch;
	name2score_type_[ "rama_prepro" ] = rama_prepro;
	name2score_type_[ "cart_bonded" ] = cart_bonded;
	name2score_type_[ "cart_bonded_angle" ] = cart_bonded_angle;
	name2score_type_[ "cart_bonded_length" ] = cart_bonded_length;
	name2score_type_[ "cart_bonded_torsion" ] = cart_bonded_torsion;
	name2score_type_[ "cart_bonded_proper" ] = cart_bonded_proper;
	name2score_type_[ "cart_bonded_improper" ] = cart_bonded_improper;
	// name2score_type_[ "csd_torsion" ] = csd_torsion; kwk commenting out csd atom type related code until I have implemented them fully
	name2score_type_[ "fa_elec" ] = fa_elec;
	name2score_type_[ "fa_elec_bb_bb" ] = fa_elec_bb_bb;
	name2score_type_[ "fa_elec_bb_sc" ] = fa_elec_bb_sc;
	name2score_type_[ "fa_elec_sc_sc" ] = fa_elec_sc_sc;
	name2score_type_[ "fa_intra_elec" ] = fa_intra_elec;
	name2score_type_[ "fa_elec_rna_phos_phos" ] = fa_elec_rna_phos_phos;
	name2score_type_[ "fa_elec_rna_phos_sugr" ] = fa_elec_rna_phos_sugr;
	name2score_type_[ "fa_elec_rna_phos_base" ] = fa_elec_rna_phos_base;
	name2score_type_[ "fa_elec_rna_sugr_sugr" ] = fa_elec_rna_sugr_sugr;
	name2score_type_[ "fa_elec_rna_sugr_base" ] = fa_elec_rna_sugr_base;
	name2score_type_[ "fa_elec_rna_base_base" ] = fa_elec_rna_base_base;

	name2score_type_[ "fa_elec_rna_phos_phos_fast" ] = fa_elec_rna_phos_phos_fast;
	name2score_type_[ "fa_elec_rna_phos_sugr_fast" ] = fa_elec_rna_phos_sugr_fast;
	name2score_type_[ "fa_elec_rna_phos_base_fast" ] = fa_elec_rna_phos_base_fast;
	name2score_type_[ "fa_elec_rna_sugr_sugr_fast" ] = fa_elec_rna_sugr_sugr_fast;
	name2score_type_[ "fa_elec_rna_sugr_base_fast" ] = fa_elec_rna_sugr_base_fast;
	name2score_type_[ "fa_elec_rna_base_base_fast" ] = fa_elec_rna_base_base_fast;

	name2score_type_[ "fa_elec_aro_aro" ] = fa_elec_aro_aro;
	name2score_type_[ "fa_elec_aro_all" ] = fa_elec_aro_all;
	name2score_type_[ "hack_aro" ] = hack_aro;
	name2score_type_[ "h2o_hbond" ] = h2o_hbond;
	name2score_type_[ "dna_dr" ] = dna_dr;
	name2score_type_[ "dna_bs" ] = dna_bs;
	name2score_type_[ "dna_bp" ] = dna_bp;
	name2score_type_[ "dna_ref" ] = dna_ref;
	name2score_type_[ "dna_env" ] = dna_env;
	name2score_type_[ "dna_pair" ] = dna_pair;
	name2score_type_[ "pro_close" ] = pro_close;
	name2score_type_[ "vdw" ] = vdw;
	name2score_type_[ "cen_hb" ] = cen_hb;
	name2score_type_[ "cenpack" ] = cenpack;
	name2score_type_[ "hybrid_vdw" ] = hybrid_vdw;
	name2score_type_[ "fa_cust_pair_dist" ] = fa_cust_pair_dist;
	name2score_type_[ "gauss" ] = gauss;
	name2score_type_[ "goap" ] = goap;
	name2score_type_[ "goap_dist" ] = goap_dist;
	name2score_type_[ "goap_angle" ] = goap_angle;

	// PyRosetta score types
#ifdef PYROSETTA
		fill_score_range(name2score_type_, "PyRosettaTwoBodyContextIndepenedentEnergy", PyRosettaTwoBodyContextIndepenedentEnergy_first, PyRosettaTwoBodyContextIndepenedentEnergy_last);
		fill_score_range(name2score_type_, "PyRosettaTwoBodyContextDependentEnergy", PyRosettaTwoBodyContextDependentEnergy_first, PyRosettaTwoBodyContextDependentEnergy_last);
		fill_score_range(name2score_type_, "PyRosettaEnergy", PyRosettaEnergy_first, PyRosettaEnergy_last);
#endif

	name2score_type_[ "python" ] = python;

	name2score_type_[ "fastsaxs" ] = fastsaxs;
	name2score_type_[ "saxs_score" ] = saxs_score;
	name2score_type_[ "saxs_fa_score" ] = saxs_fa_score;
	name2score_type_[ "saxs_cen_score" ] = saxs_cen_score;
	name2score_type_[ "pddf_score" ] = pddf_score;

	name2score_type_[ "fiberdiffraction" ] = fiberdiffraction;
	name2score_type_[ "fiberdiffractiondens" ] = fiberdiffractiondens;
#ifdef USECUDA
	name2score_type_[ "fiberdiffractiongpu" ] = fiberdiffractiongpu;
#endif

	name2score_type_[ "fa_pair" ] = fa_pair; // fa_pair == fa_pair_pol_pol
	name2score_type_[ "fa_pair_aro_aro" ] = fa_pair_aro_aro;
	name2score_type_[ "fa_pair_aro_pol" ] = fa_pair_aro_pol;
	name2score_type_[ "fa_pair_pol_pol" ] = fa_pair_pol_pol;
	name2score_type_[ "fa_plane" ] = fa_plane;
	name2score_type_[ "hbond_sr_bb" ] = hbond_sr_bb;
	name2score_type_[ "hbond_lr_bb" ] = hbond_lr_bb;
	name2score_type_[ "hbond_bb_sc" ] = hbond_bb_sc;
	name2score_type_[ "hbond_sr_bb_sc" ] = hbond_sr_bb_sc;
	name2score_type_[ "hbond_lr_bb_sc" ] = hbond_lr_bb_sc;
	name2score_type_[ "hbond_sc"    ] = hbond_sc;
	name2score_type_[ "hbond"    ] = hbond;
	name2score_type_[ "fa_grpelec" ] = fa_grpelec;
	name2score_type_[ "interchain_pair" ] = interchain_pair;
	name2score_type_[ "interchain_vdw" ] = interchain_vdw;
	name2score_type_[ "interface_dd_pair" ] = interface_dd_pair;

	name2score_type_[ "ch_bond"    ] = ch_bond;
	name2score_type_[ "ch_bond_bb_bb" ] = ch_bond_bb_bb;
	name2score_type_[ "ch_bond_sc_sc" ] = ch_bond_sc_sc;
	name2score_type_[ "ch_bond_bb_sc" ] = ch_bond_bb_sc;

	name2score_type_[ "neigh_vect"  ] = neigh_vect;
	name2score_type_[ "neigh_count" ] = neigh_count;
	name2score_type_[ "neigh_vect_raw"] = neigh_vect_raw;
	name2score_type_[ "symE_bonus"  ] = symE_bonus;
	name2score_type_[ "sym_lig"  ] = sym_lig;

	name2score_type_["orbitals_hpol_bb"] = orbitals_hpol_bb;
	name2score_type_["pci_cation_pi"] = pci_cation_pi;
	name2score_type_["pci_pi_pi"] =pci_pi_pi;
	name2score_type_["pci_salt_bridge"] =pci_salt_bridge;
	name2score_type_["pci_hbond"] =pci_hbond;


	name2score_type_[ "geom_sol"    ] = geom_sol;
	name2score_type_[ "occ_sol_fitted"    ] = occ_sol_fitted;
	name2score_type_[ "occ_sol_fitted_onebody"    ] = occ_sol_fitted_onebody;
	name2score_type_[ "occ_sol_exact"    ] = occ_sol_exact;

	name2score_type_[ "gb_elec" ] = gb_elec;
	name2score_type_[ "multipole_elec" ] = multipole_elec;
	name2score_type_[ "fa_sasa" ] = fa_sasa;

	name2score_type_[ "facts_elec" ] = facts_elec;
	name2score_type_[ "facts_solv" ] = facts_solv;
	name2score_type_[ "facts_sasa" ] = facts_sasa;

	name2score_type_[ "PB_elec" ] = PB_elec;
	name2score_type_[ "dslf_ss_dst" ] = dslf_ss_dst;
	name2score_type_[ "dslf_cs_ang" ] = dslf_cs_ang;
	name2score_type_[ "dslf_ss_dih" ] = dslf_ss_dih;
	name2score_type_[ "dslf_ca_dih" ] = dslf_ca_dih;
	name2score_type_[ "dslf_cbs_ds" ] = dslf_cbs_ds;
	name2score_type_[ "dslf_fa13" ] = dslf_fa13;
	name2score_type_[ "dslfc_cen_dst" ] = dslfc_cen_dst;
	name2score_type_[ "dslfc_cb_dst"  ] = dslfc_cb_dst;
	name2score_type_[ "dslfc_ang"     ] = dslfc_ang;
	name2score_type_[ "dslfc_cb_dih"  ] = dslfc_cb_dih;
	name2score_type_[ "dslfc_bb_dih"  ] = dslfc_bb_dih;

	name2score_type_[ "dslfc_rot"  ] = dslfc_rot;
	name2score_type_[ "dslfc_trans"  ] = dslfc_trans;
	name2score_type_[ "dslfc_RT"  ] = dslfc_RT;

	name2score_type_[ "custom_atom_pair" ] = custom_atom_pair;
	name2score_type_[ "atom_pair_constraint" ] = atom_pair_constraint;
	name2score_type_[ "base_pair_constraint" ] = base_pair_constraint;
	name2score_type_[ "coarse_chainbreak_constraint" ] = coarse_chainbreak_constraint;
	name2score_type_[ "dunbrack_constraint" ] = dunbrack_constraint;
	name2score_type_[ "angle_constraint" ] = angle_constraint;
	name2score_type_[ "dihedral_constraint" ] = dihedral_constraint;
	name2score_type_[ "big_bin_constraint" ] = big_bin_constraint;
	name2score_type_[ "constant_constraint" ] = constant_constraint;
	name2score_type_[ "coordinate_constraint" ] = coordinate_constraint;
	name2score_type_[ "site_constraint" ] = site_constraint;
	name2score_type_[ "metalhash_constraint" ] = metalhash_constraint;
	name2score_type_[ "metalbinding_constraint" ] = metalbinding_constraint;
	name2score_type_[ "bond_geometry"] = bond_geometry;
	name2score_type_[ "Hpol_bond_geometry"] = Hpol_bond_geometry;

	name2score_type_[ "rama"    ] = rama;
	name2score_type_[ "rama2b"  ] = rama2b;
	name2score_type_[ "omega"    ] = omega;
	name2score_type_[ "fa_dun" ] = fa_dun;
	name2score_type_[ "fa_dun_dev" ] = fa_dun_dev;
	name2score_type_[ "fa_dun_rot" ] = fa_dun_rot;
	name2score_type_[ "fa_dun_semi" ] = fa_dun_semi;
	name2score_type_[ "dna_chi" ] = dna_chi;
	name2score_type_[ "p_aa_pp" ] = p_aa_pp;
	name2score_type_[ "p_aa_ss" ] = p_aa_ss;
	name2score_type_[ "yhh_planarity" ] = yhh_planarity;
	name2score_type_[ "hxl_tors" ] = hxl_tors;
	name2score_type_[ "h2o_intra" ] =  h2o_intra;
	name2score_type_[ "ref" ] = ref;
	name2score_type_[ "ref_nc" ] = ref_nc;
	name2score_type_[ "seqdep_ref" ] = seqdep_ref;
	name2score_type_[ "nmer_ref" ] = nmer_ref;
	name2score_type_[ "nmer_pssm" ] = nmer_pssm;
	name2score_type_[ "nmer_svm" ] = nmer_svm;
	name2score_type_[ "envsmooth" ] = envsmooth;
	name2score_type_[ "e_pH" ] = e_pH;
	name2score_type_[ "rna_bulge"] = rna_bulge;

	// Context-Independent, One-Body Carbohydrate Scoring Terms
	name2score_type_[ "sugar_bb"] = sugar_bb;

	name2score_type_[ "loop_close"] = loop_close;
	name2score_type_[ "missing_res"] = missing_res;
	name2score_type_[ "free_suite"] = free_suite;
	name2score_type_[ "free_2HOprime"] = free_2HOprime;
	name2score_type_[ "free_side_chain"] = free_side_chain;
	name2score_type_[ "free_base"] = free_base;
	name2score_type_[ "free_res"] = free_res;
	name2score_type_[ "free_dof"] = free_dof;
	name2score_type_[ "intermol"] = intermol;
	// Variant type to flag rotamers for alternative scoring with varying weight
	name2score_type_[ "special_rot"] = special_rot;
	name2score_type_[ "dna_dihedral_bb"] = dna_dihedral_bb;
	name2score_type_[ "dna_dihedral_chi"] = dna_dihedral_chi;
	name2score_type_[ "dna_dihedral_sugar"] = dna_dihedral_sugar;

	name2score_type_[ "other_pose"] = other_pose;

	name2score_type_[ "env" ]    = env;
	name2score_type_[ "burial" ] = burial; //Probably unused
	name2score_type_[ "burial_v2" ] = burial_v2; //burial score from a set of residues
	name2score_type_[ "abego" ]  = abego;
	name2score_type_[ "pair" ]   = pair;
	name2score_type_[ "cbeta" ]  = cbeta;
	name2score_type_[ "DFIRE" ]  = DFIRE;

	//Centroid Rotamer Model
	name2score_type_[ "cen_rot_pair" ] = cen_rot_pair;
	name2score_type_[ "cen_rot_pair_ang" ] = cen_rot_pair_ang;
	name2score_type_[ "cen_rot_pair_dih" ] = cen_rot_pair_dih;
	name2score_type_[ "cen_rot_env" ] = cen_rot_env;
	name2score_type_[ "cen_rot_dun" ] = cen_rot_dun;
	name2score_type_[ "cen_rot_cbeta" ] = cen_rot_cbeta;

	//bw membrane scoring terms
	name2score_type_[ "Menv" ] = Menv;
	name2score_type_[ "Menv_non_helix" ] = Menv_non_helix;
	name2score_type_[ "Menv_termini" ] = Menv_termini;
	name2score_type_[ "Menv_tm_proj" ] = Menv_tm_proj;
	name2score_type_[ "Mcbeta" ] = Mcbeta;
	name2score_type_[ "Mpair" ] = Mpair;
	name2score_type_[ "Mlipo" ] = Mlipo;
	//pba membrane all atom terms
	name2score_type_[ "fa_mbenv" ] = fa_mbenv;
	name2score_type_[ "fa_mbsolv" ] = fa_mbsolv;
	name2score_type_[ "Menv_smooth" ] = Menv_smooth;

	// Membrane Framework SupportedEnergy Terms
	// centroid - added by @ralford 3/30/14
	name2score_type_[ "mp_env" ] = MPEnv;
	name2score_type_[ "mp_cbeta" ] = MPCbeta;
	name2score_type_[ "mp_pair" ] = MPPair;
	name2score_type_[ "mp_lipo" ] = MPLipo;
	name2score_type_[ "mp_termini" ] = MPTermini;
	name2score_type_[ "mp_nonhelix" ] = MPNonHelix;
	name2score_type_[ "mp_tmproj" ] = MPTMProj;

	// fullatom - added by @ralford 5/14/14
	name2score_type_[ "fa_mpenv" ] = FaMPEnv;
	name2score_type_[ "fa_mpsolv" ] = FaMPSolv;
	name2score_type_[ "fa_mpenv_smooth" ] = FaMPEnvSmooth;

	name2score_type_[ "rg" ] = rg;
	name2score_type_[ "rg_local" ] = rg_local;
	name2score_type_[ "co" ] = co;
	name2score_type_[ "peptide_bond" ] = peptide_bond;
	name2score_type_[ "pcs" ] = pcs;
	name2score_type_["pcsTs1"]=pcsTs1;
	name2score_type_["pcsTs2"]=pcsTs2;
	name2score_type_["pcsTs3"]=pcsTs3;
	name2score_type_["pcsTs4"]=pcsTs4;
	name2score_type_[ "pcs2" ] = pcs2;
	name2score_type_[ "dock_ens_conf" ] = dock_ens_conf;

	//fpd smooth (differentiable) centroid terms
	name2score_type_[ "cen_hb" ] = cen_hb;
	name2score_type_[ "cen_env_smooth" ] = cen_env_smooth;
	name2score_type_[ "cen_pair_smooth" ] = cen_pair_smooth;
	name2score_type_[ "cbeta_smooth" ] = cbeta_smooth;
	name2score_type_[ "cenpack_smooth" ] = cenpack_smooth;
	name2score_type_[ "hs_pair" ] = hs_pair;
	name2score_type_[ "ss_pair" ] = ss_pair;
	name2score_type_[ "rsigma" ] = rsigma;
	name2score_type_[ "sheet" ] = sheet;
	name2score_type_[ "csa" ] = csa;
	name2score_type_[ "dc" ] = dc;
	name2score_type_[ "rdc" ] = rdc;
	name2score_type_[ "rdc_segments" ] = rdc_segments;
	name2score_type_[ "rdc_rohl" ] =rdc_rohl;
	name2score_type_[ "holes" ] = holes;
	name2score_type_[ "holes_resl" ] = holes_resl;
	name2score_type_[ "holes_decoy" ] = holes_decoy;
	name2score_type_[ "holes_min" ] = holes_min;
	name2score_type_[ "holes_min_mean" ] = holes_min_mean;
	name2score_type_[ "rna_chem_shift" ] = rna_chem_shift;
	name2score_type_[ "rna_chem_map"] = rna_chem_map;
	name2score_type_[ "rna_chem_map_lores"] = rna_chem_map_lores;
	name2score_type_[ "dab_sasa" ] = dab_sasa;
	name2score_type_[ "dab_sev" ] = dab_sev;
	name2score_type_[ "sa" ] = sa;
	name2score_type_[ "d2h_sa" ] = d2h_sa;
	name2score_type_[ "ProQM" ] = ProQM;
	name2score_type_[ "ProQ" ] = ProQ;

	//Centroid motif derived score functions
	name2score_type_[ "cen_pair_motifs"] = cen_pair_motifs;
	name2score_type_[ "cen_pair_motif_degree"] = cen_pair_motif_degree;
	name2score_type_[ "ss_contact_worst"] = ss_contact_worst;

	name2score_type_[ "interchain_env"] = interchain_env;
	name2score_type_[ "interchain_contact"] = interchain_contact;

	name2score_type_[ "rna_rg"] = rna_rg;
	name2score_type_[ "rna_vdw"] = rna_vdw;
	name2score_type_[ "rna_base_backbone"] = rna_base_backbone;
	name2score_type_[ "rna_backbone_backbone"] = rna_backbone_backbone;
	name2score_type_[ "rna_repulsive"] = rna_repulsive;

	name2score_type_[ "rna_base_pair"] = rna_base_pair;
	name2score_type_[ "rna_base_axis"] = rna_base_axis;
	name2score_type_[ "rna_base_stagger"] = rna_base_stagger;
	name2score_type_[ "rna_base_stack"] = rna_base_stack;
	name2score_type_[ "rna_base_stack_axis"] = rna_base_stack_axis;

	name2score_type_[ "rna_data_base"] = rna_data_base;
	name2score_type_[ "rna_data_backbone"] = rna_data_backbone;

	name2score_type_[ "rna_mg_point"] = rna_mg_point;
	name2score_type_[ "rna_mg_point_indirect"] = rna_mg_point_indirect;

	name2score_type_[ "mg"] = mg;
	name2score_type_[ "mg_lig"] = mg_lig;
	name2score_type_[ "mg_sol"] = mg_sol;
	name2score_type_[ "mg_ref"] = mg_ref;
	name2score_type_[ "hoh_ref"] = hoh_ref;

	//Will these ever really be used?
	name2score_type_[ "rna_base_pair_pairwise"] = rna_base_pair_pairwise;
	name2score_type_[ "rna_base_axis_pairwise"] = rna_base_axis_pairwise;
	name2score_type_[ "rna_base_stagger_pairwise"] = rna_base_stagger_pairwise;
	name2score_type_[ "rna_base_stack_pairwise"] = rna_base_stack_pairwise;
	name2score_type_[ "rna_base_stack_axis_pairwise"] = rna_base_stack_axis_pairwise;

	name2score_type_[ "fa_stack"]       = fa_stack;
	name2score_type_[ "fa_stack_lower"] = fa_stack_lower;
	name2score_type_[ "fa_stack_upper"] = fa_stack_upper;
	name2score_type_[ "fa_stack_aro"]   = fa_stack_aro;
	name2score_type_[ "fa_stack_ext"]   = fa_stack_ext;
	name2score_type_[ "fa_stack_sol"]   = fa_stack_sol;
	name2score_type_[ "fa_stack_lr"]    = fa_stack_lr;

	name2score_type_[ "stack_elec"] = stack_elec;
	name2score_type_[ "stack_elec_base_base"] = stack_elec_base_base;
	name2score_type_[ "stack_elec_base_bb"] = stack_elec_base_bb;

	name2score_type_[ "rna_torsion"] = rna_torsion;
	name2score_type_[ "rna_torsion_sc"] = rna_torsion_sc;
	name2score_type_[ "rna_suite"] = rna_suite;
	name2score_type_[ "rna_jr_suite"] = rna_jr_suite;
	name2score_type_[ "suiteness_bonus"] = suiteness_bonus;
	name2score_type_[ "rna_sugar_close"] = rna_sugar_close;
	name2score_type_[ "rna_bond_geometry"] = rna_bond_geometry;
	name2score_type_[ "rna_fa_atr_base"] = rna_fa_atr_base;
	name2score_type_[ "rna_fa_rep_base"] = rna_fa_rep_base;

	////////////Intra-res RNA specific score terms//////////////////
	name2score_type_[ "lk_polar_intra_RNA" ] = lk_polar_intra_RNA;
	name2score_type_[ "lk_nonpolar_intra_RNA" ] = lk_nonpolar_intra_RNA;
	name2score_type_[ "fa_intra_RNA_base_phos_atr" ] = fa_intra_RNA_base_phos_atr;
	name2score_type_[ "fa_intra_RNA_base_phos_rep" ] = fa_intra_RNA_base_phos_rep;
	name2score_type_[ "fa_intra_RNA_base_phos_sol" ] = fa_intra_RNA_base_phos_sol;
	name2score_type_[ "hbond_intra" ] = hbond_intra; //Currently affects only RNA.
	name2score_type_[ "geom_sol_intra_RNA"    ] = geom_sol_intra_RNA;
	name2score_type_[ "geom_sol_fast"    ] = geom_sol_fast;
	name2score_type_[ "geom_sol_fast_intra_RNA"    ] = geom_sol_fast_intra_RNA;

	name2score_type_[ "dna_bb_torsion"] = dna_bb_torsion;
	name2score_type_[ "dna_sugar_close"] = dna_sugar_close;
	name2score_type_[ "dna_base_distance"] = dna_base_distance;

	name2score_type_[ "chainbreak" ] = chainbreak;
	name2score_type_[ "linear_chainbreak" ] = linear_chainbreak;
	name2score_type_[ "overlap_chainbreak" ] = overlap_chainbreak;
	name2score_type_[ "distance_chainbreak" ] = distance_chainbreak;
	name2score_type_[ "dof_constraint" ] = dof_constraint;
	name2score_type_[ "rms_energy" ] = rms;
	name2score_type_[ "suck" ] = suck;
	name2score_type_[ "res_type_constraint" ] = res_type_constraint;
	name2score_type_[ "res_type_linking_constraint" ] = res_type_linking_constraint;
	name2score_type_[ "pocket_constraint" ] = pocket_constraint;
	name2score_type_[ "backbone_stub_constraint" ] = backbone_stub_constraint;
	name2score_type_[ "backbone_stub_linear_constraint" ] = backbone_stub_linear_constraint;
	name2score_type_[ "pack_stat" ] = pack_stat;

	name2score_type_[ "surface" ] = surface;
	name2score_type_[ "hpatch" ] = hpatch;
	name2score_type_[ "p_aa" ] = p_aa;
	name2score_type_[ "unfolded" ] = unfolded;
	name2score_type_[ "split_unfolded_two_body" ] = split_unfolded_two_body;

	// individual terms for split unfolded weight calculation
	name2score_type_[ "fa_atr_ref" ] = fa_atr_ref;
	name2score_type_[ "fa_rep_ref" ] = fa_rep_ref;
	name2score_type_[ "fa_sol_ref" ] = fa_sol_ref;
	name2score_type_[ "fa_elec_ref" ] = fa_elec_ref;
	name2score_type_[ "hbond_ref" ] = hbond_ref;
	name2score_type_[ "dslf_fa13_ref" ] = dslf_fa13_ref;
	name2score_type_[ "fa_intra_atr_ref" ] = fa_intra_atr_ref;
	name2score_type_[ "fa_intra_rep_ref" ] = fa_intra_rep_ref;
	name2score_type_[ "fa_intra_sol_ref" ] = fa_intra_sol_ref;
	name2score_type_[ "pro_close_ref" ] = pro_close_ref;
	name2score_type_[ "fa_dun_ref" ] = fa_dun_ref;
	name2score_type_[ "fa_dun_dev_ref" ] = fa_dun_dev_ref;
	name2score_type_[ "fa_dun_rot_ref" ] = fa_dun_rot_ref;
	name2score_type_[ "fa_dun_semi_ref" ] = fa_dun_semi_ref;
	name2score_type_[ "rama_ref" ] = rama_ref;
	name2score_type_[ "p_aa_pp_ref" ] = p_aa_pp_ref;
	name2score_type_[ "omega_ref" ] = omega_ref;
	name2score_type_[ "mm_lj_intra_rep_ref" ] = mm_lj_intra_rep_ref;
	name2score_type_[ "mm_lj_intra_atr_ref" ] = mm_lj_intra_atr_ref;
	name2score_type_[ "mm_twist_ref" ] = mm_twist_ref;

	name2score_type_[ "elec_dens_fast" ] = elec_dens_fast;
	name2score_type_[ "elec_dens_window" ] = elec_dens_window;
	name2score_type_[ "elec_dens_whole_structure_ca" ] = elec_dens_whole_structure_ca;
	name2score_type_[ "elec_dens_whole_structure_allatom" ] = elec_dens_whole_structure_allatom;
	name2score_type_[ "elec_dens_atomwise" ] = elec_dens_atomwise;
	name2score_type_[ "grid_vdw" ] = grid_vdw;
	name2score_type_[ "xtal_ml" ] = xtal_ml;
	name2score_type_[ "xtal_rwork" ] = xtal_rwork;
	name2score_type_[ "xtal_rfree" ] = xtal_rfree;

	name2score_type_[ "natbias_ss" ] = natbias_ss;
	name2score_type_[ "natbias_hs" ] = natbias_hs;
	name2score_type_[ "natbias_hh" ] = natbias_hh;
	name2score_type_[ "natbias_stwist" ] = natbias_stwist;

	name2score_type_[ "aa_cmp" ] = aa_cmp;

	name2score_type_[ "ring_close" ] = ring_close; //General score term for enforcing ring closure in proline-like noncanonicals.  NOTE: EITHER ring_close, OR pro_close, OR cart_bonded should be used -- otherwise we'll double-count!
	name2score_type_[ "aa_repeat" ] = aa_repeat; //A wholebody score term for penalizing long stretches of repeat sequence (e.g. poly-Q sequences).
	name2score_type_[ "aa_composition" ] = aa_composition; //A wholebody score term for penalizing deviation from a desired amino acid composition.
	name2score_type_[ "aspartimide_penalty"] = aspartimide_penalty; //A context-independent two-body score term for penalizing two-residue sequences likely to produce the aspartimide side-product during peptide synthesis.

	name2score_type_[ "sidechain_neighbors" ] = sidechain_neighbors;

	name2score_type_[ "total_score" ] = total_score;

	name2score_type_[ "dummy_score_type" ] = dummy_score_type;

	runtime_assert( name2score_type_.size() == end_of_score_type_enumeration );

	score_type2name_.resize( end_of_score_type_enumeration );
	for ( std::map< std::string, ScoreType >::const_iterator iter = name2score_type_.begin(),
			iter_end = name2score_type_.end(); iter != iter_end; ++iter ) {
		score_type2name_[ iter->second ] = iter->first;
	}

}


//////////////////////////////////////////////////////////////////////////////
/// @brief give a ScoreType string name and return its enum type
ScoreType
ScoreTypeManager::score_type_from_name( std::string const & name )
{
	setup_score_type_names();
	std::map< std::string, ScoreType >::const_iterator iter( name2score_type_.find( name ) );
	if ( iter == name2score_type_.end() ) {
		utility_exit_with_message("unrecognized score_type type "+name);
	}
	return iter->second;
}

std::string
ScoreTypeManager::name_from_score_type( ScoreType score_type )
{
	setup_score_type_names();
	return score_type2name_[ score_type ];
}

/// @brief
bool
ScoreTypeManager::is_score_type( std::string const & name )
{
	setup_score_type_names();
	std::map< std::string, ScoreType >::const_iterator iter( name2score_type_.find( name ) );
	return iter != name2score_type_.end();
}

}
}
