// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/FullatomDisulfidePotential.cc
/// @brief  Fullatom Disulfide potential class definition
/// @author Bill Schief
/// @author blindly ported by Andrew Leaver-Fay
/// @author rewritten by Robert Vernon

// Unit Headers
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>

// Package Headers
#include <core/scoring/disulfides/DisulfideAtomIndices.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <basic/database/open.hh>

// Utility Headers
#include <numeric/xyz.functions.hh>
#include <utility/exit.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/interpolation/Histogram.hh>
#include <numeric/statistics.functions.hh>

#include <numeric/deriv/angle_deriv.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.disulfides.FullatomDisulfidePotential");

namespace core {
namespace scoring {
namespace disulfides {

using namespace core;

//-----------------------------------------
// physical parameters of disulfide bonds
//-----------------------------------------
//
// (Computed from disulfide containing proteins in the 2008 Vall database)
//
// CB--SG--SG angle (gaussian)
//      CBSGSG MEAN:  1.819120   SDEV: 0.104496
//
// CB--SG--SG--CB dihedral angle (bimodal gaussian)
// posCBSGSGCB MEAN:  1.641427   SDEV: 0.247741
// negCBSGSGCB MEAN: -1.517302   SDEV: 0.203992
//
// CA--CB--SG--SG dihedral angles (coarse bimodal gaussian)
// posCACBSGSG MEAN:  1.709486   SDEV: 0.660339
// negCACBSGSG MEAN: -1.471579   SDEV: 0.604909
//
// Notes:
// CBSGSG fits well to a gaussian curve
//
// CBSGSGCB is sharply bimodal across a ~180 degree angle difference
// with a negative angle peak that's ~24% larger than the positive angle peak.
// We may actually want to upweight the larger peak, but the 'true difference' between
// the peaks correlates with a ton of other angles and energies so for now they can compete equally
//
// CACBSGSG is coarsely bimodal and not a true gaussian.
// There are peaks ~140 degrees removed from each other on either side of zero,
// but the two regions between the peaks differ in behavior
// (specifically: spanning zero is forbidden, spanning -180/180 is allowed)
// For now the score ignores that and treats them as two gaussians with large standard deviations.
//
//
// These angles are specific to the following functions:
//
//	CACBSGSG Dihedral = cbsg_dihedral_func_( new CBSG_Dihedral_Func ) = dslf_ss_dih
//	CBSGSGCB Dihedral = sgsg_dihedral_func_( new SGSG_Dihedral_Func ) = dslf_ca_dih
//	    CBSGSG Aangle =      cb_angle_func_( new CB_Angle_Func )      = dslf_cs_ang
//
//  SG--SG Distances are done by the following rosetta++ antique, which uses a fitted histogram
//
//  sg_dist_func_( new SG_Dist_Func )
//
//  -Robert Vernon (March 2010)



//
//	core::Real const ideal_ss_dist_in_disulfide = { 2.02 };
//	 // mean sulfur-sulfur distance in natives // angstroms
//
//	core::Real const disulf_ss_dist_stdev = { 0.35 };
//	 // standard dev. of s-s dist in natives // degrees 0.35
//
//	core::Real const ideal_cs_angle_in_disulfide = { 103.4 };
//	 // mean cbeta-sulfur-sulfur angle in natives // degrees
//
//	core::Real const disulf_cs_angle_stdev = { 5.0 };
//	 // standard dev. of cbeta-s-s angle in natives // degrees  2.6
//
//	core::Real const mean_pos_dihedral_in_disulf = { 87.9 };
//	 // mean positive cbeta-sulfur-sulfur-cbeta dihedral angle // degrees
//
//	core::Real const disulf_pos_dihedral_stdev = { 21.8 };
//	 // standard dev. of pos. dihedral angle in natives // degrees
//
//	core::Real const mean_pos_dihedral_gauss1 = { 87.4 };
//	 // mean positive cbeta-sulfur-sulfur-cbeta dihedral angle // degrees
//
//	core::Real const stdev_pos_dihedral_gauss1 = { 20.9 };
//	 // standard dev. of pos. dihedral angle in natives // degrees
//
//	core::Real const mean_pos_dihedral_gauss2 = { 95.6 };
//	 // mean positive cbeta-sulfur-sulfur-cbeta dihedral angle // degrees
//
//	core::Real const stdev_pos_dihedral_gauss2 = { 3.0 };
//	 // standard dev. of pos. dihedral angle in natives // degrees
//
//	core::Real const mean_neg_dihedral_in_disulf = { -86.2 };
//	 // mean negative cbeta-sulfur-sulfur-cbeta dihedral angle // degrees
//
//	core::Real const disulf_neg_dihedral_stdev = { 11.1 };
//	 // standard dev. of neg. dihedral angle in natives // degrees
//
//	core::Real const ideal_ca_dihedral_in_disulf = { 74.0 };
//	 // ideal calpha-cbeta-sulfur-sulfur dihedral (abs val) // degrees


FullatomDisulfidePotential::FullatomDisulfidePotential() :
	parent(),
	disulf_ssdist_cutoff_( 3.0 ),
	cbsg_dihedral_func_( new CBSG_Dihedral_Func ),
	sgsg_dihedral_func_( new SGSG_Dihedral_Func ),
	cb_angle_func_( new CB_Angle_Func ),
	sg_dist_func_( new SG_Dist_Func ),
	wt_dihSS_(0.1),
	wt_dihCS_(0.1),
	wt_ang_(0.1),
	wt_len_(0.1),
	shift_(2.0)
{
	mest_ = exp(-20.0);
}

FullatomDisulfidePotential::~FullatomDisulfidePotential() {}

void
FullatomDisulfidePotential::print_score_functions() const
{
	using namespace numeric::constants::d;
	for (Real angle = -360.0; angle <= 360.0; angle = angle + 0.1) {

		std::cout << "ANGLE: " << angle
							<< " CBSGSG " << cb_angle_func_->func(angle*degrees_to_radians)
							<< " CBSGSGCB " << sgsg_dihedral_func_->func(angle*degrees_to_radians)
							<< " CACBSGSG " << cbsg_dihedral_func_->func(angle*degrees_to_radians)
							<< std::endl;
	}

	for (Real distance = 0.0; distance <= 10; distance = distance + 0.01) {
		std::cout << "DISTANCE: " << distance << " SCORE: " << sg_dist_func_->func(distance) << std::endl;
	}
}

/**
 * @brief Calculated scores for a disulfide bond between two residues
 *
 * @param[in] res1 The lower residue of the disulfide
 * @param[in] res2 The upper residue of the disulfide. Assumed to be bonded to res1
 * @param[out] distance_score_this_disulfide A score based on S-S distance
 * @param[out] csangles_score_this_disulfide A score based on the Cb-S angles
 * @param[out] dihedral_score_this_disulfide A score based on the S-S dihedral
 * @param[out] ca_dihedral_sc_this_disulf A score based on the Cb-S dihedrals
 * @param[out] truefalse_fa_disulf True if these residues should be disulfide bonded
 *
 * @details Scores are interpolated from the histograms created by the
 *  farlx_*_initializer methods.
 *  The distance score has such a histogram as its core with two linear
 *  functions at either side so that the score increases to infinity
 */
void
FullatomDisulfidePotential::score_this_disulfide_old(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	DisulfideAtomIndices const & res1_atom_indices,
	DisulfideAtomIndices const & res2_atom_indices,
	Energy & distance_score_this_disulfide,
	Energy & csangles_score_this_disulfide,
	Energy & dihedral_score_this_disulfide,
	Energy & ca_dihedral_sc_this_disulf,
	bool & truefalse_fa_disulf
	) const
{
	// Allocate memory for disulf params. Values are set by get_disulfide_params
	// dist between cys sulfurs
	Real ssdist(-1);
	// Cb-S-S bond angles
	Real csang_1(-1),csang_2(-1);
	// dihedral (torsion) angle, cbeta1-s1-s2-cbeta2
	Real dihed(360);
	// dihedral (torsion) angle, calpha1-cbeta1-s1-s2
	Real disulf_ca_dihedral_angle_1(360);
	// dihedral (torsion) angle, calpha2-cbeta2-s2-s1
	Real disulf_ca_dihedral_angle_2(360);

	//Assume res1 and res2 are bonded and score them
	get_disulfide_params(res1,res2,res1_atom_indices,res2_atom_indices,
		ssdist,csang_1,csang_2,dihed,disulf_ca_dihedral_angle_1,disulf_ca_dihedral_angle_2);

	distance_score_this_disulfide = 0;
	csangles_score_this_disulfide = 0;
	dihedral_score_this_disulfide = 0;
	ca_dihedral_sc_this_disulf = 0;
	truefalse_fa_disulf = false;

	// ssdist score:

	distance_score_this_disulfide = sg_dist_func_->func(ssdist);

	// csangle score:

	//Check that csang is in the right domain
	runtime_assert_msg( ! ( csang_1 > 180. || csang_2 > 180. ) , "Error csang > 180" );

	using namespace numeric::constants::d;
	Energy csang_1_score = cb_angle_func_->func(csang_1*degrees_to_radians);
	Energy csang_2_score = cb_angle_func_->func(csang_2*degrees_to_radians);
	csangles_score_this_disulfide = (csang_1_score+csang_2_score)*0.5;

	// dihedral score:

	dihedral_score_this_disulfide = sgsg_dihedral_func_->func(dihed*degrees_to_radians);

	// ca_dihedral score: (this is just the chi_2 angle, rotation about cbeta-sulfur)

	// chi_2 of cys_1

	ca_dihedral_sc_this_disulf = cbsg_dihedral_func_->func(disulf_ca_dihedral_angle_1*degrees_to_radians);
	ca_dihedral_sc_this_disulf += cbsg_dihedral_func_->func(disulf_ca_dihedral_angle_2*degrees_to_radians);
	//Compensate for adding the scores together
	ca_dihedral_sc_this_disulf *= .5 ;


	//std::cout << "DUMP_CA_DIH " << disulf_ca_dihedral_angle_1 << " " << cbsg_dihedral_func_->func(disulf_ca_dihedral_angle_1*degrees_to_radians) << std::endl;
	//std::cout << "DUMP_CA_DIH " << disulf_ca_dihedral_angle_2 << " " << cbsg_dihedral_func_->func(disulf_ca_dihedral_angle_2*degrees_to_radians) << std::endl;

	//std::cout << "DUMP_SS_DIH " << dihed << " " << sgsg_dihedral_func_->func(dihed*degrees_to_radians) << std::endl;

	//std::cout << "DUMP_ANG " << csang_1 << " " << cb_angle_func_->func(csang_1*degrees_to_radians) << std::endl;
	//	std::cout << "DUMP_ANG " << csang_2 << " " << cb_angle_func_->func(csang_2*degrees_to_radians) << std::endl;

	//std::cout << "DUMP_DIST " << ssdist << " " << sg_dist_func_->func(ssdist) << std::endl;

	// Call it a disulfide or not!
	if ( ssdist < disulf_ssdist_cutoff_ ) {
		truefalse_fa_disulf = true;
	}
}


void
FullatomDisulfidePotential::get_disulfide_derivatives_old(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	DisulfideAtomIndices const & res1_atom_indices,
	DisulfideAtomIndices const & res2_atom_indices,
	Size const at1,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace id;
	using namespace constraints;

	ResiduePairXYZ respairxyz( res1, res2 );
	Vector f1( 0.0 ), f2( 0.0 );
	if ( res1_atom_indices.derivative_atom( at1 ) == CYS_C_ALPHA ) {

		f1 = f2 = 0;
		DihedralConstraint dihedral_ang_cst = DihedralConstraint(
			AtomID( at1, res1.seqpos() ),
			AtomID( res1_atom_indices.c_beta_index(), res1.seqpos() ),
			AtomID( res1_atom_indices.disulf_atom_index(), res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
				cbsg_dihedral_func_, dslf_ca_dih );
		dihedral_ang_cst.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

	} else if ( res1_atom_indices.derivative_atom( at1 ) == CYS_C_BETA  ) {

		f1 = f2 = 0;
		AngleConstraint ang_cst = AngleConstraint(
			AtomID( at1, res1.seqpos() ),
			AtomID( res1_atom_indices.disulf_atom_index(), res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
				cb_angle_func_, dslf_cs_ang );
		ang_cst.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		DihedralConstraint dihedral_ang_cst = DihedralConstraint(
			AtomID( res1_atom_indices.c_alpha_index(), res1.seqpos() ),
			AtomID( at1, res1.seqpos() ),
			AtomID( res1_atom_indices.disulf_atom_index(), res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
				cbsg_dihedral_func_, dslf_ca_dih );
		dihedral_ang_cst.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		DihedralConstraint ss_dihedral_ang_cst = DihedralConstraint(
			AtomID( at1, res1.seqpos() ),
			AtomID( res1_atom_indices.disulf_atom_index(), res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
			AtomID( res2_atom_indices.c_beta_index(), res2.seqpos() ),
				sgsg_dihedral_func_, dslf_ss_dih );
		ss_dihedral_ang_cst.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;

	} else if ( res1_atom_indices.derivative_atom( at1 ) == CYS_S_GAMMA ) {

		AtomPairConstraint apc( AtomID( at1, res1.seqpos() ), AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ), sg_dist_func_, dslf_ss_dst );
		apc.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;

		f1 = f2 = 0;
		AngleConstraint ang_cst1 = AngleConstraint(
			AtomID( res1_atom_indices.c_beta_index(), res1.seqpos() ),
			AtomID( at1, res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
				cb_angle_func_, dslf_cs_ang );
		ang_cst1.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		AngleConstraint ang_cst2 = AngleConstraint(
			AtomID( at1, res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
			AtomID( res2_atom_indices.c_beta_index(), res2.seqpos() ),
				cb_angle_func_, dslf_cs_ang );
		ang_cst2.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		DihedralConstraint dihedral_ang_cst1 = DihedralConstraint(
			AtomID( res1_atom_indices.c_alpha_index(), res1.seqpos() ),
			AtomID( res1_atom_indices.c_beta_index(), res1.seqpos() ),
			AtomID( at1, res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
				cbsg_dihedral_func_, dslf_ca_dih );
		dihedral_ang_cst1.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		DihedralConstraint dihedral_ang_cst2 = DihedralConstraint(
			AtomID( res2.atom_index("CA"), res2.seqpos() ),
			AtomID( res2_atom_indices.c_beta_index(), res2.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
			AtomID( at1, res1.seqpos() ),
				cbsg_dihedral_func_, dslf_ca_dih );
		dihedral_ang_cst2.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += 0.5 * f1;
		F2 += 0.5 * f2;

		f1 = f2 = 0;
		DihedralConstraint ss_dihedral_ang_cst = DihedralConstraint(
			AtomID( res1_atom_indices.c_beta_index(), res1.seqpos() ),
			AtomID( at1, res1.seqpos() ),
			AtomID( res2_atom_indices.disulf_atom_index(), res2.seqpos() ),
			AtomID( res2_atom_indices.c_beta_index(), res2.seqpos() ),
				sgsg_dihedral_func_, dslf_ss_dih );
		ss_dihedral_ang_cst.fill_f1_f2( AtomID( at1, res1.seqpos() ), respairxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;

	}
}


void
FullatomDisulfidePotential::score_this_disulfide(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	DisulfideAtomIndices const & res1_atom_indices,
	DisulfideAtomIndices const & res2_atom_indices,
	Energy & score
) const {
	using numeric::constants::f::pi;
	using namespace numeric::statistics;
	using namespace core::chemical;

	Real ssdist, csang_1, csang_2, dihed, disulf_ca_dihedral_angle_1, disulf_ca_dihedral_angle_2;
	get_disulfide_params(res1,res2,res1_atom_indices,res2_atom_indices,
		ssdist,csang_1,csang_2,dihed,disulf_ca_dihedral_angle_1,disulf_ca_dihedral_angle_2);

	score = -shift_;

	{ // distance
		// z <- (x-location)/scale;
	  // score <- x^2/2 - Log[Erfc[-((s x)/Sqrt[2])]] + (1/2) (Log[2] + Log[\[Pi]]) + Log[s]
		core::Real z = (ssdist-params_.d_location)/params_.d_scale;
		core::Real score_d = z*z/2 - log( errfc( -params_.d_shape*z / sqrt(2.) ) + mest_ );
		score += wt_len_*score_d;
	}

	{ // angles
		Real ang_d = csang_1;
		Real score_a = -params_.a_logA - params_.a_kappa*cos( pi/180 * (ang_d-params_.a_mu) );
		score += wt_ang_*score_a;

		ang_d = csang_2;
		score_a = -params_.a_logA - params_.a_kappa*cos( pi/180 * (ang_d-params_.a_mu) );
		score += wt_ang_*score_a;
	}

	{ // SS dih
		Real ang_ss = dihed;
		Real exp_score1 = exp(params_.dss_logA1)*exp(params_.dss_kappa1*cos(  pi/180 * (ang_ss-params_.dss_mu1) ));
		Real exp_score2 = exp(params_.dss_logA2)*exp(params_.dss_kappa2*cos(  pi/180 * (ang_ss-params_.dss_mu2) ));
		Real score_ss = -log (exp_score1 + exp_score2 + mest_);
		score += wt_dihSS_*score_ss;
	}

	{ // CB-S dihedrals
		Real ang_cs = disulf_ca_dihedral_angle_1;
		Real exp_score1 = exp(params_.dcs_logA1)*exp(params_.dcs_kappa1*cos(  pi/180 * (ang_cs-params_.dcs_mu1) ));
		Real exp_score2 = exp(params_.dcs_logA2)*exp(params_.dcs_kappa2*cos(  pi/180 * (ang_cs-params_.dcs_mu2) ));
		Real exp_score3 = exp(params_.dcs_logA3)*exp(params_.dcs_kappa3*cos(  pi/180 * (ang_cs-params_.dcs_mu3) ));
		Real score_cs = -log (exp_score1 + exp_score2 + exp_score3 + mest_);
		score += wt_dihCS_*score_cs;

		ang_cs = disulf_ca_dihedral_angle_2;
		exp_score1 = exp(params_.dcs_logA1)*exp(params_.dcs_kappa1*cos(  pi/180 * (ang_cs-params_.dcs_mu1) ));
		exp_score2 = exp(params_.dcs_logA2)*exp(params_.dcs_kappa2*cos(  pi/180 * (ang_cs-params_.dcs_mu2) ));
		exp_score3 = exp(params_.dcs_logA3)*exp(params_.dcs_kappa3*cos(  pi/180 * (ang_cs-params_.dcs_mu3) ));
		score_cs = -log (exp_score1 + exp_score2 + exp_score3 + mest_);
		score += wt_dihCS_*score_cs;
	}
}


void
FullatomDisulfidePotential::get_disulfide_derivatives(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	DisulfideAtomIndices const & res1_atom_indices,
	DisulfideAtomIndices const & res2_atom_indices,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {
	using numeric::constants::f::pi;
	using namespace numeric::statistics;
	using namespace core::chemical;

	Real ssdist, csang_1, csang_2, dihed, disulf_ca_dihedral_angle_1, disulf_ca_dihedral_angle_2;
	get_disulfide_params(res1,res2,res1_atom_indices,res2_atom_indices,
		ssdist,csang_1,csang_2,dihed,disulf_ca_dihedral_angle_1,disulf_ca_dihedral_angle_2);

	// atom indices
	Size i_ca1=res1_atom_indices.c_alpha_index(), i_cb1=res1_atom_indices.c_beta_index(), i_sg1=res1_atom_indices.disulf_atom_index();
	Size i_ca2=res2_atom_indices.c_alpha_index(), i_cb2=res2_atom_indices.c_beta_index(), i_sg2=res2_atom_indices.disulf_atom_index();

	// storage
	Vector f1,f2;
	Real d, theta, phi;

 	{ // distance
		// z <- (x-location)/scale;
	  // score <- x^2/2 - Log[Erfc[-((s x)/Sqrt[2])]] + (1/2) (Log[2] + Log[\[Pi]]) + Log[s]
		core::Real z = (ssdist-params_.d_location)/params_.d_scale;
		core::Real dscore_d = z/params_.d_scale -
			( exp( -0.5*z*z*params_.d_shape*params_.d_shape ) * sqrt(2./pi) * params_.d_shape ) / (params_.d_scale * errfc(-params_.d_shape*z / sqrt(2.) ) + 1.e-12 );
		dscore_d = weights[ dslf_fa13 ]*wt_len_*dscore_d;

		numeric::deriv::distance_f1_f2_deriv( res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), d, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_d * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_d * f2;
		r2_atom_derivs[ i_sg2 ].f1() -= dscore_d * f1;
		r2_atom_derivs[ i_sg2 ].f2() -= dscore_d * f2;
	}

	{ // angles
		Real ang_d = csang_1;
		Real dscore_a = params_.a_kappa * sin( pi/180 * (ang_d-params_.a_mu) );
		dscore_a = weights[ dslf_fa13 ]*wt_ang_*dscore_a;
		numeric::deriv::angle_p1_deriv( res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), theta, f1, f2 );
		r1_atom_derivs[ i_cb1 ].f1() += dscore_a * f1;
		r1_atom_derivs[ i_cb1 ].f2() += dscore_a * f2;
		numeric::deriv::angle_p2_deriv( res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), theta, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_a * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_a * f2;
		numeric::deriv::angle_p1_deriv( res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), res1.xyz( i_cb1 ), theta, f1, f2 );
		r2_atom_derivs[ i_sg2 ].f1() += dscore_a * f1;
		r2_atom_derivs[ i_sg2 ].f2() += dscore_a * f2;

		ang_d = csang_2;
		dscore_a = params_.a_kappa * sin( pi/180 * (ang_d-params_.a_mu) );
		dscore_a = weights[ dslf_fa13 ]*wt_ang_*dscore_a;
		numeric::deriv::angle_p1_deriv( res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), theta, f1, f2 );
		r2_atom_derivs[ i_cb2 ].f1() += dscore_a * f1;
		r2_atom_derivs[ i_cb2 ].f2() += dscore_a * f2;
		numeric::deriv::angle_p2_deriv( res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), theta, f1, f2 );
		r2_atom_derivs[ i_sg2 ].f1() += dscore_a * f1;
		r2_atom_derivs[ i_sg2 ].f2() += dscore_a * f2;
		numeric::deriv::angle_p1_deriv( res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), res2.xyz( i_cb2 ), theta, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_a * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_a * f2;
	}

	{ // SS dih
		Real ang_ss = dihed;
		Real exp_score1 = exp(params_.dss_logA1)*exp(params_.dss_kappa1*cos(  pi/180 * (ang_ss-params_.dss_mu1) ));
		Real exp_score2 = exp(params_.dss_logA2)*exp(params_.dss_kappa2*cos(  pi/180 * (ang_ss-params_.dss_mu2) ));
		Real dscore_ss = 0.0;
		dscore_ss += exp_score1 * params_.dss_kappa1 * sin( pi/180 * (ang_ss-params_.dss_mu1) );
		dscore_ss += exp_score2 * params_.dss_kappa2 * sin( pi/180 * (ang_ss-params_.dss_mu2) );
		dscore_ss /= ( exp_score1 + exp_score2 + mest_);
		dscore_ss = weights[ dslf_fa13 ]*wt_dihSS_*dscore_ss;

		numeric::deriv::dihedral_p1_cosine_deriv( res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), res2.xyz( i_cb2 ), phi, f1, f2 );
		r1_atom_derivs[ i_cb1 ].f1() += dscore_ss * f1;
		r1_atom_derivs[ i_cb1 ].f2() += dscore_ss * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), res2.xyz( i_cb2 ), phi, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_ss * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_ss * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), res1.xyz( i_cb1 ), phi, f1, f2 );
		r2_atom_derivs[ i_sg2 ].f1() += dscore_ss * f1;
		r2_atom_derivs[ i_sg2 ].f2() += dscore_ss * f2;
		numeric::deriv::dihedral_p1_cosine_deriv( res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), res1.xyz( i_cb1 ), phi, f1, f2 );
		r2_atom_derivs[ i_cb2 ].f1() += dscore_ss * f1;
		r2_atom_derivs[ i_cb2 ].f2() += dscore_ss * f2;
	}

	{ // CB-S dihedrals
		Real ang_cs = disulf_ca_dihedral_angle_1;
		Real exp_score1 = exp(params_.dcs_logA1)*exp(params_.dcs_kappa1*cos(  pi/180 * (ang_cs-params_.dcs_mu1) ));
		Real exp_score2 = exp(params_.dcs_logA2)*exp(params_.dcs_kappa2*cos(  pi/180 * (ang_cs-params_.dcs_mu2) ));
		Real exp_score3 = exp(params_.dcs_logA3)*exp(params_.dcs_kappa3*cos(  pi/180 * (ang_cs-params_.dcs_mu3) ));
		Real dscore_cs = 0.0;
		dscore_cs += exp_score1 * params_.dcs_kappa1 * sin( pi/180 * (ang_cs-params_.dcs_mu1) );
		dscore_cs += exp_score2 * params_.dcs_kappa2 * sin( pi/180 * (ang_cs-params_.dcs_mu2) );
		dscore_cs += exp_score3 * params_.dcs_kappa3 * sin( pi/180 * (ang_cs-params_.dcs_mu3) );
		dscore_cs /= ( exp_score1 + exp_score2 + exp_score3 + mest_);
		dscore_cs = weights[ dslf_fa13 ]*wt_dihCS_*dscore_cs;

		numeric::deriv::dihedral_p1_cosine_deriv( res1.xyz( i_ca1 ), res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), phi, f1, f2 );
		r1_atom_derivs[ i_ca1 ].f1() += dscore_cs * f1;
		r1_atom_derivs[ i_ca1 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res1.xyz( i_ca1 ), res1.xyz( i_cb1 ), res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), phi, f1, f2 );
		r1_atom_derivs[ i_cb1 ].f1() += dscore_cs * f1;
		r1_atom_derivs[ i_cb1 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), res1.xyz( i_cb1 ), res1.xyz( i_ca1 ), phi, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_cs * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p1_cosine_deriv( res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), res1.xyz( i_cb1 ), res1.xyz( i_ca1 ), phi, f1, f2 );
		r2_atom_derivs[ i_sg2 ].f1() += dscore_cs * f1;
		r2_atom_derivs[ i_sg2 ].f2() += dscore_cs * f2;


		ang_cs = disulf_ca_dihedral_angle_2;
		exp_score1 = exp(params_.dcs_logA1)*exp(params_.dcs_kappa1*cos(  pi/180 * (ang_cs-params_.dcs_mu1) ));
		exp_score2 = exp(params_.dcs_logA2)*exp(params_.dcs_kappa2*cos(  pi/180 * (ang_cs-params_.dcs_mu2) ));
		exp_score3 = exp(params_.dcs_logA3)*exp(params_.dcs_kappa3*cos(  pi/180 * (ang_cs-params_.dcs_mu3) ));
		dscore_cs = 0.0;
		dscore_cs += exp_score1 * params_.dcs_kappa1 * sin( pi/180 * (ang_cs-params_.dcs_mu1) );
		dscore_cs += exp_score2 * params_.dcs_kappa2 * sin( pi/180 * (ang_cs-params_.dcs_mu2) );
		dscore_cs += exp_score3 * params_.dcs_kappa3 * sin( pi/180 * (ang_cs-params_.dcs_mu3) );
		dscore_cs /= ( exp_score1 + exp_score2 + exp_score3 + mest_);
		dscore_cs = weights[ dslf_fa13 ]*wt_dihCS_*dscore_cs;

		numeric::deriv::dihedral_p1_cosine_deriv( res2.xyz( i_ca2 ), res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), phi, f1, f2 );
		r2_atom_derivs[ i_ca2 ].f1() += dscore_cs * f1;
		r2_atom_derivs[ i_ca2 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res2.xyz( i_ca2 ), res2.xyz( i_cb2 ), res2.xyz( i_sg2 ), res1.xyz( i_sg1 ), phi, f1, f2 );
		r2_atom_derivs[ i_cb2 ].f1() += dscore_cs * f1;
		r2_atom_derivs[ i_cb2 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p2_cosine_deriv( res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), res2.xyz( i_cb2 ), res2.xyz( i_ca2 ), phi, f1, f2 );
		r2_atom_derivs[ i_sg2 ].f1() += dscore_cs * f1;
		r2_atom_derivs[ i_sg2 ].f2() += dscore_cs * f2;
		numeric::deriv::dihedral_p1_cosine_deriv( res1.xyz( i_sg1 ), res2.xyz( i_sg2 ), res2.xyz( i_cb2 ), res2.xyz( i_ca2 ), phi, f1, f2 );
		r1_atom_derivs[ i_sg1 ].f1() += dscore_cs * f1;
		r1_atom_derivs[ i_sg1 ].f2() += dscore_cs * f2;
	}
}


///////////////////////////////////////////
/// Private: Methods, Data Initializers ///
///////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// @brief Calculates several geometrical parameters of a fullatom disulfide bond
/// @details
///   given residue num for 2 cys involved in disulfide bond, returns
///   four quantities:
///   sulf-sulf dist, 2 carb-sulf bond angles, and a 4-atom dihedral angle.
///   Angles are returned in degrees.
///
/// @param[in]   coord1 - in - fullatom coords of cys 1
/// @param[in]   coord2 - in - fullatom coords of cys 2
/// @param[out]   dist_between_sulfurs - out - s1-s2 distance
/// @param[out]   cs_bond_angle_1 - out - cb1-s1-s2 bond angle
/// @param[out]   cs_bond_angle_2 - out - cb2-s1-s2 bond angle
/// @param[out]   disulf_dihedral_angle - out - cb1-s1-s2-cb2 dihedral
/// @param[out]   disulf_ca_dihedral_angle_1 - out - ca1-cb1-s1-s2 dihedral
/// @param[out]   disulf_ca_dihedral_angle_2 - out - ca2-cb2-s2-s1 dihedral
///
/// @author Bill Schief
/////////////////////////////////////////////////////////////////////////////////
void
FullatomDisulfidePotential::get_disulfide_params(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	DisulfideAtomIndices const & res1_atom_indices,
	DisulfideAtomIndices const & res2_atom_indices,
	Distance & dist_between_sulfurs, // dist between cys sulfurs
	Real & cs_bond_angle_1,
	Real & cs_bond_angle_2,
	Real & disulf_dihedral_angle, // dihedral (torsion) angle, cbeta-s-s-cbeta
	Real & disulf_ca_dihedral_angle_1,
	 // dihedral (torsion) angle, calpha1-cbeta1-s1-s2
	Real & disulf_ca_dihedral_angle_2 // dihedral (torsion) angle, calpha2-cbeta2-s2-s1
) const
{
	using namespace numeric::constants::d;

	Vector calpha_1( res1.atom( res1_atom_indices.c_alpha_index() ).xyz() );
	Vector cbeta_1 ( res1.atom( res1_atom_indices.c_beta_index( ) ).xyz() );
	Vector sulfur_1( res1.atom( res1_atom_indices.disulf_atom_index() ).xyz() );
	Vector calpha_2( res2.atom( res2_atom_indices.c_alpha_index() ).xyz() );
	Vector cbeta_2 ( res2.atom( res2_atom_indices.c_beta_index( ) ).xyz() );
	Vector sulfur_2( res2.atom( res2_atom_indices.disulf_atom_index() ).xyz() );

	dist_between_sulfurs = sulfur_1.distance( sulfur_2 );
	cs_bond_angle_1 = angle_of( cbeta_1, sulfur_1, sulfur_2);
	cs_bond_angle_2 = angle_of( cbeta_2, sulfur_2, sulfur_1);
	cs_bond_angle_1 *= radians_to_degrees; // convert
	cs_bond_angle_2 *= radians_to_degrees; // convert
	disulf_dihedral_angle      = dihedral_degrees(cbeta_1,sulfur_1,sulfur_2,cbeta_2);
	disulf_ca_dihedral_angle_1 = dihedral_degrees(calpha_1,cbeta_1,sulfur_1,sulfur_2);
	disulf_ca_dihedral_angle_2 = dihedral_degrees(calpha_2,cbeta_2,sulfur_2,sulfur_1);
}


//------------------------------------------------------------------------------

CBSG_Dihedral_Func::CBSG_Dihedral_Func() :
	//cbsg_pos_peak_( 86.0 ),
	//cbsg_pos_sd_( 72.0 ),
	//cbsg_neg_peak_( -74.0 ),
	//cbsg_neg_sd_( 39.0 )

	//csf_cbang1a_(0.1,0.9,10,0.0),
	//csf_cbang1b_(0.1,0.9,10,0.0)

	csf_cbang1_(0.1,0.9,10,0.0),
	csf_cbang2_(1.35,1.1,-5,0.0),
	csf_cbang3_(-1.25,1.5,-5,0.0)
{}

CBSG_Dihedral_Func::~CBSG_Dihedral_Func() {}

/// @param ang[in] Dihedral angle in radians
/// @note This function is not continuous across 0 or 180 degrees. This is bad,
///  and is only acceptable because it occurs at the maximums so we minimize
///  away from the discontinuity. -Spencer
Real
CBSG_Dihedral_Func::func( Real const ang ) const
{
// 	using namespace numeric::constants::d;

// 	Real ang_deg = ang * radians_to_degrees;
// 	Real const ang_peak  = ang_deg > 0.0 ? cbsg_pos_peak_ : cbsg_neg_peak_;
// 	Real const ang_sd    = ang_deg > 0.0 ? cbsg_pos_sd_   : cbsg_neg_sd_;
// 	Real const delta_ang = ang_deg - ang_peak;
// 	Real const ang_frac  = delta_ang / ang_sd;
// 	return -std::exp(-(ang_frac*ang_frac));
	//return csf_cbang1a_.func(ang) + csf_cbang1b_.func(ang);

	return csf_cbang1_.func(ang) + (csf_cbang2_.func(ang))/10 + (csf_cbang3_.func(ang))/10 + 0.095;
}

/// @param ang[in] Dihedral angle in radians
Real
CBSG_Dihedral_Func::dfunc( Real const ang ) const {
// 	using namespace numeric::constants::d;

// 	Real ang_deg = ang * radians_to_degrees;
// 	Real const ang_peak  = ang_deg > 0.0 ? cbsg_pos_peak_ : cbsg_neg_peak_;
// 	Real const ang_sd    = ang_deg > 0.0 ? cbsg_pos_sd_   : cbsg_neg_sd_;
// 	Real const delta_ang = ang_deg - ang_peak;
// 	Real const ang_frac  = delta_ang / ang_sd;

// 	return radians_to_degrees * std::exp(-(ang_frac*ang_frac)) * ( 2*ang_deg - 2*ang_peak ) / ( ang_sd * ang_sd );
	return  csf_cbang1_.dfunc(ang) + (csf_cbang2_.dfunc(ang))/10 + (csf_cbang3_.dfunc(ang))/10;
}

///////////////

// CB--SG--SG--CB dihedral angle (bimodal gaussian)
// posCBSGSGCB MEAN:  1.641427   SDEV: 0.247741
// negCBSGSGCB MEAN: -1.517302   SDEV: 0.203992

SGSG_Dihedral_Func::SGSG_Dihedral_Func() :
	//csf_cbang1a_( 1.641426, 0.25, -2, 0.0),
	//csf_cbang2a_( 1.641426, 0.7, -40, 0.0),
	//csf_cbang1b_(-1.517302, 0.25, -2, 0.0),
	//csf_cbang2b_(-1.517302, 0.7, -40, 0.0)
csf_cbang1a_( 1.641426, 0.25, -2.3, 0.0),
csf_cbang2a_( 1.641426, 0.9, -20, 0.0),
csf_cbang1b_(-1.517302, 0.25, -2.3, 0.0),
csf_cbang2b_(-1.517302, 0.9, -20, 0.0)
 {}

SGSG_Dihedral_Func::~SGSG_Dihedral_Func() {}

/// @param ang[in] Dihedral angle in radians
Real
SGSG_Dihedral_Func::func( Real const ang) const
{
// 	using namespace numeric::constants::d;

// 	Real ang_deg = ang * radians_to_degrees;
// 	Real score(0);
// 	fa_sgsg_dihedral_scores()->interpolate(ang_deg,score);
// 	return score;

	//std::cout << "SGSG " << ang << " " << csf_cbang1a_.func(ang) << " " << csf_cbang2a_.func(ang) << " " << csf_cbang1b_.func(ang) << " " << csf_cbang2b_.func(ang) << " " << 2*csf_cbang1a_.func(ang) + csf_cbang2a_.func(ang) + 2*csf_cbang1b_.func(ang) + csf_cbang2b_.func(ang) << std::endl;

	//return 2*csf_cbang1a_.func(ang) + csf_cbang2a_.func(ang) + 2*csf_cbang1b_.func(ang) + csf_cbang2b_.func(ang) + 0.6602761;

	//return 10*csf_cbang1a_.func(ang) + csf_cbang2a_.func(ang) + 10*csf_cbang1b_.func(ang) + csf_cbang2b_.func(ang);
	return 25*csf_cbang1a_.func(ang) + csf_cbang2a_.func(ang) + 25*csf_cbang1b_.func(ang) + csf_cbang2b_.func(ang) + 4.58;
}

/// @param ang[in] Dihedral angle in radians
Real
SGSG_Dihedral_Func::dfunc( Real const ang ) const
{
// 	using namespace numeric::constants::d;

// 	Real ang_deg = ang * radians_to_degrees;
// 	Real deriv;
// 	fa_sgsg_dihedral_scores()->derivative(ang_deg, deriv); // units/degree
// 	return deriv*radians_to_degrees; // return units/radian

//	return 10*csf_cbang1a_.dfunc(ang) + csf_cbang2a_.dfunc(ang) + 10*csf_cbang1b_.dfunc(ang) + csf_cbang2b_.dfunc(ang);
	return 25*csf_cbang1a_.dfunc(ang) + csf_cbang2a_.dfunc(ang) + 25*csf_cbang1b_.dfunc(ang) + csf_cbang2b_.dfunc(ang);
}

/// Access the histogram for this Func (in degrees)
// numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
// SGSG_Dihedral_Func::fa_sgsg_dihedral_scores()
// {
// 	using namespace numeric::interpolation;
// 	static HistogramCOP<Real,Real>::Type scores(0);
// 	if(scores == 0) {
// 		utility::io::izstream scores_stream;
// 		basic::database::open( scores_stream, "scoring/score_functions/disulfides/fa_CbSSCb_dihedral_score");
// 		scores = new Histogram<Real,Real>( scores_stream() );
// 		scores_stream.close();
// 	}
// 	return scores;
// }

///////////////

CB_Angle_Func::CB_Angle_Func() :
	//chf_cbang_( 1.819120, 0.104496 )
	//chf_cbang_( 1.819120, 0.4 ),

	//onno
	//csf_cbang_( 1.819120, 0.1, -10, 0.0)
	//ommo

	//CHECKIN
	//csf_cbang_( 1.819120, 0.104496, -200, 0.0)

	csf_cbang1_( 1.819120, 0.208961, -200, 0.0),
	csf_cbang2_( 1.819120, 0.208961, -50, 0.0)

//csf_cbang1_( 1.819120, 0.3, -2, 0.0),
	//csf_cbang2_( 1.819120, 0.1, -10, 0.0)
	//2c
													//csf_cbang1_( 1.819120, 0.104496, -7, 0.0),
														//csf_cbang2_( 1.819120, 0.104496, -0.5, 0.0)
	//csf_cbang1_( 1.819120, 0.104496, -200 ),
	//csf_cbang2_( 1.819120, 0.2, -20 ),
	//csf_cbang3_( 1.819120, 0.2, -5, 0.0)
{}

CB_Angle_Func::~CB_Angle_Func() {}

/// @param csang[in] Dihedral angle in radians
Real
CB_Angle_Func::func( Real const ang) const
{
// 	using namespace numeric::constants::d;

// 	Real ang_deg = ang*radians_to_degrees;
// 	//Check that csang is in the right domain
// 	runtime_assert_msg(  0 <= ang_deg && ang_deg <= 180., "Error csang > 180" );

// 	Real score(0);
// 	// Will be out of the histogram's range sometimes, but we like the default
// 	// behavior of using the boundary scores for extreme ang_deg
// 	CB_Angle_Func::fa_csang_scores()->interpolate(ang_deg,score);
// 	return score;

	//std::cout << "CBANG " << ang << " " << csf_cbang1_.func(ang) << " " << csf_cbang2_.func(ang) << " " << csf_cbang3_.func(ang) << " " << (0.2*chf_cbang_.func(ang)) + csf_cbang1_.func(ang) + csf_cbang2_.func(ang) + csf_cbang3_.func(ang) << std::endl;

	//return (0.2*chf_cbang_.func(ang)) + csf_cbang1_.func(ang) + csf_cbang2_.func(ang) + csf_cbang3_.func(ang);

	//std::cout << "CBANG " << ang << " " << csf_cbang1_.func(ang) << " " << csf_cbang2_.func(ang) << " " << (50*csf_cbang1_.func(ang)) + (1000.0*csf_cbang2_.func(ang)) << std::endl;


	//return (csf_cbang1_.func(ang)) + (20.0*csf_cbang2_.func(ang)) - 0.393803696;
	//return (csf_cbang1_.func(ang)) + (csf_cbang2_.func(ang)) - 0.393803696;
	//return (csf_cbang1_.func(ang)) + (csf_cbang2_.func(ang)) - 0.393803696;

	//onno
	//return 10*(csf_cbang_.func(ang)) + std::abs(ang - 1.819120) + 2.44918662;

	//CHECKIN
	//return 2*(csf_cbang_.func(ang)) + std::abs(ang - 1.819120) + 2.0;

	return (csf_cbang1_.func(ang)) + (csf_cbang2_.func(ang)) + 4*(std::abs(ang - 1.819120)) + 2.0;
}

/// @param ang[in] Dihedral angle in radians
Real
CB_Angle_Func::dfunc( Real const ang ) const
{
// 	using namespace numeric::constants::d;
// 	Real ang_deg = ang * radians_to_degrees;
// 	//std::cout << "ang_deg: " << ang_deg << std::endl;

// 	Real d_csang_score_dang( 0.0 );
// 	// Should be zero when ang_deg is not in the hist range
// 	CB_Angle_Func::fa_csang_scores()->derivative(ang_deg,d_csang_score_dang);
// 	return d_csang_score_dang * radians_to_degrees;
//	return (0.2*chf_cbang_.dfunc(ang)) + csf_cbang1_.dfunc(ang) + csf_cbang2_.dfunc(ang) + csf_cbang3_.dfunc(ang);

	//return (csf_cbang1_.dfunc(ang)) + (20.0*csf_cbang2_.dfunc(ang));

	//onno
	//return 10*(csf_cbang_.dfunc(ang)) + (ang - 1.819120)/std::abs(ang - 1.819120);

	//CHECKIN
	//return 2*(csf_cbang_.dfunc(ang)) + (ang - 1.819120)/std::abs(ang - 1.819120);

	return (csf_cbang1_.dfunc(ang)) + (csf_cbang2_.dfunc(ang)) + 4*((ang - 1.819120)/std::abs(ang - 1.819120));
}

/// Access the histogram for this Func (in degrees)
// numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
// CB_Angle_Func::fa_csang_scores()
// {
// 	using namespace numeric::interpolation;
// 	static HistogramCOP<Real,Real>::Type scores(0);
// 	if(scores == 0) {
// 		utility::io::izstream scores_stream;
// 		basic::database::open( scores_stream, "scoring/score_functions/disulfides/fa_CaCbS_angle_score");
// 		scores = new Histogram<Real,Real>( scores_stream() );
// 		scores_stream.close();
// 	}
// 	return scores;
// }

//////////////////

SG_Dist_Func::SG_Dist_Func() {}

SG_Dist_Func::~SG_Dist_Func() {}

/// @param ssdist[in] S-S distance, in Angstroms
Real
SG_Dist_Func::func( Real const ssdist ) const
{
	Real distance_score_this_disulfide( 0.0 );

	numeric::interpolation::HistogramCOP<Real,Real>::Type scores =
		SG_Dist_Func::fa_ssdist_scores();
	Real ssdist_min_range( scores->minimum() );
	Real ssdist_max_range( scores->maximum() );

	if ( ssdist >= ssdist_max_range ) { // max tail
		Real ssdist_score_at_max(0);
		scores->interpolate(ssdist_max_range,ssdist_score_at_max);
		Real ssdist_max_tail_slope(0);
		scores->derivative(ssdist_max_range-scores->step_size(),ssdist_max_tail_slope);

		distance_score_this_disulfide = ssdist_score_at_max + ( ssdist - ssdist_max_range ) * ssdist_max_tail_slope;
	} else if ( ssdist <= ssdist_min_range ) { // min tail
		Real ssdist_score_at_min(0);
		scores->interpolate(ssdist_min_range,ssdist_score_at_min);
		Real ssdist_min_tail_slope(0);
		scores->derivative(ssdist_min_range,ssdist_min_tail_slope);

		distance_score_this_disulfide = ssdist_score_at_min + ( ssdist - ssdist_min_range ) * ssdist_min_tail_slope;
	} else {
		SG_Dist_Func::fa_ssdist_scores()->interpolate(ssdist, distance_score_this_disulfide);
	}

	return distance_score_this_disulfide;
}

/// @param ssdist[in] S-S distance, in Angstroms
Real
SG_Dist_Func::dfunc( Real const ssdist ) const {
	Real d_distance_score_this_disulfide_ddis( 0.0 );

	numeric::interpolation::HistogramCOP<Real,Real>::Type scores =
		SG_Dist_Func::fa_ssdist_scores();
	Real ssdist_min_range( scores->minimum() );
	Real ssdist_max_range( scores->maximum() );
	if ( ssdist >= ssdist_max_range ) { // max tail
		scores->derivative(ssdist_max_range-scores->step_size(),d_distance_score_this_disulfide_ddis);
	} else if ( ssdist <= ssdist_min_range ) { // min tail
		scores->derivative(ssdist_min_range,d_distance_score_this_disulfide_ddis);
	} else {
		SG_Dist_Func::fa_ssdist_scores()->derivative(ssdist, d_distance_score_this_disulfide_ddis);
	}
	return d_distance_score_this_disulfide_ddis;
}

/// Access the histogram for this Func
numeric::interpolation::HistogramCOP<core::Real,core::Real>::Type
SG_Dist_Func::fa_ssdist_scores()
{
	using namespace numeric::interpolation;
	static HistogramCOP<Real,Real>::Type scores(0);
	if(scores == 0) {
		utility::io::izstream scores_stream;
		basic::database::open( scores_stream, "scoring/score_functions/disulfides/fa_SS_distance_score");
		scores = new Histogram<Real,Real>( scores_stream() );
		scores_stream.close();
	}
	return scores;
}

}
}
}

