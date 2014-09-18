// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/CentroidDisulfidePotential.cc
/// @brief  Centroid Disulfide Energy Potentials
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date   12/17/08

// Unit Headers
#include <core/scoring/disulfides/CentroidDisulfidePotential.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/database/open.hh>
#include <core/scoring/constraints/util.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/constants.hh>
#include <numeric/interpolation/Histogram.hh>


#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>

static thread_local basic::Tracer TR( "core.scoring.disulfides.CentroidDisulfidePotential" );


using namespace core;
using core::scoring::disulfides::CentroidDisulfidePotential;
using namespace numeric::interpolation;
using core::conformation::Residue;
using std::string;
using utility::vector1;

namespace core {
namespace scoring {
namespace disulfides {

/**
 * Constructor
 */
CentroidDisulfidePotential::CentroidDisulfidePotential() {}

/**
 * Deconstructor
 */
CentroidDisulfidePotential::~CentroidDisulfidePotential() {}

/**
 * @brief Calculates scoring terms for the disulfide bond specified
 * @note Equivalent to the expanded form, but discards geometry info
 */
void
CentroidDisulfidePotential::score_disulfide(
		Residue const & res1,
		Residue const & res2,
		Energy & cbcb_distance_score,
		Energy & centroid_distance_score,
		Energy & cacbcb_angle_1_score,
		Energy & cacbcb_angle_2_score,
		Energy & cacbcbca_dihedral_score,
		Energy & backbone_dihedral_score
		) const
{
	Real cbcb_distance_sq, centroid_distance_sq, cacbcb_angle_1, cacbcb_angle_2,
	cacbcbca_dihedral, backbone_dihedral, score_factor;
	score_disulfide(res1,
			res2,
			cbcb_distance_sq,
			centroid_distance_sq,
			cacbcb_angle_1,
			cacbcb_angle_2,
			cacbcbca_dihedral,
			backbone_dihedral,
			cbcb_distance_score,
			centroid_distance_score,
			cacbcb_angle_1_score,
			cacbcb_angle_2_score,
			cacbcbca_dihedral_score,
			backbone_dihedral_score,
			score_factor
			);
}

/**
 * @brief Calculates scoring terms and geometry
 *
 * If a full atom pose is given, centroid_distance_score will be zero.
 *
 * If one of the residues is glycine it will be replaced with alanine for
 * the scores which require a CB atom.
 *
 * @note distances are given squared to avoid calling std::sqrt unnecessarily
 */
void
CentroidDisulfidePotential::score_disulfide(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real & cbcb_distance_sq,
		core::Real & centroid_distance_sq,
		core::Real & cacbcb_angle_1,
		core::Real & cacbcb_angle_2,
		core::Real & cacbcbca_dihedral,
		core::Real & backbone_dihedral,
		core::Energy & cbcb_distance_score,
		core::Energy & centroid_distance_score,
		core::Energy & cacbcb_angle_1_score,
		core::Energy & cacbcb_angle_2_score,
		core::Energy & cacbcbca_dihedral_score,
		core::Energy & backbone_dihedral_score,
		core::Real & cb_score_factor
		) const
{
	//The range of cb distances present in nature (squared)
	static const Real min_native_cb_dist_sq = 10; //ang^2
	static const Real max_native_cb_dist_sq = 22; //ang^2
	//Cutoff to calculate the angle terms
	static const Real max_cb_dist_sq = 400; //ang^2

	//Calculate the distances and angles of this disulfide
	disulfide_params(res1, res2,
			cbcb_distance_sq,
			centroid_distance_sq,
			cacbcb_angle_1,
			cacbcb_angle_2,
			cacbcbca_dihedral,
			backbone_dihedral);


	//Interpolate scores from the parameters
	//Do the unscaled scores here, then the reweighted scores
	cbcb_distance_score       = cb_distance_func_->func(cbcb_distance_sq);

	if(centroid_distance_sq < 0) {
		//Indicates error, probably full atom mode
		centroid_distance_score = 0.0;
	} else {
		centroid_distance_score   = cen_distance_func_->func(centroid_distance_sq);
	}

	//Score factor: Reweight angle scores based on the cb distance squared
	//0-10A^2       0.0
	//10-22A^2      1.0
	//22-400A^2     linear function between (22,1) and (400,0)
	//>400A^2       0.0
	if( cbcb_distance_sq < min_native_cb_dist_sq || max_cb_dist_sq < cbcb_distance_sq ) {
		//cb_score_factor = 0.; don't bother scoring
		cacbcb_angle_1_score      = 0.;
		cacbcb_angle_2_score      = 0.;
		cacbcbca_dihedral_score   = 0.;
		backbone_dihedral_score   = 0.;
		return;
	}
	else if( cbcb_distance_sq < max_native_cb_dist_sq ) { //native range
		cb_score_factor = 1.;
	}
	else { //longer than native, so linearly decay to zero
		//slope = 1/( 22 - 400), x-intercept = 400
		cb_score_factor = (cbcb_distance_sq - max_cb_dist_sq) /
			(max_native_cb_dist_sq - max_cb_dist_sq);
	}

	cacbcb_angle_1_score      = cacbcb_angle_func_->func(cacbcb_angle_1);
	cacbcb_angle_2_score      = cacbcb_angle_func_->func(cacbcb_angle_2);
	cacbcbca_dihedral_score   = cacbcbca_dihedral_func_->func(cacbcbca_dihedral);
	backbone_dihedral_score   = ncacac_dihedral_func_->func(backbone_dihedral);

	cacbcbca_dihedral_score *= cb_score_factor;
	backbone_dihedral_score *= cb_score_factor;
	cacbcb_angle_1_score *= cb_score_factor;
	cacbcb_angle_2_score *= cb_score_factor;
}


/**
 * @brief calculates some degrees of freedom between two centroid cys residues
 *
 * If one of the residues is glycine it will be substituted with an idealize
 * alanine geometry for the calculations which require a Cb molecule.
 *
 * centroid_distance requires CEN atoms be defined. If full atom residues
 * are specified this function returns centroid_distance of -1.
 *
 * @param cbcb_distance     The distance between Cbetas squared
 * @param centroid_distance The distance between centroids squared
 * @param cacbcb_angle_1    The Ca1-Cb1-Cb2 planar angle, in degrees
 * @param cacbcb_angle_2    The Ca2-Cb2-Cb1 planar angle, in degrees
 * @param cacbcbca_dihedral The Ca1-Cb1-Cb2-Ca2 dihedral angle
 * @param backbone_dihedral The N-Ca1-Ca2-C2 dihedral angle
 */
void
CentroidDisulfidePotential::disulfide_params(
		Residue const& res1,
		Residue const& res2,
		Real & cbcb_distance_sq,
		Real & centroid_distance_sq,
		Real & cacbcb_angle_1,
		Real & cacbcb_angle_2,
		Real & cacbcbca_dihedral,
		Real & backbone_dihedral)
{
	using numeric::constants::d::radians_to_degrees;

	conformation::ResidueCAP res1_ptr(&res1);
	conformation::ResidueCAP res2_ptr(&res2);
	// Glycines pose a problem because they have no CB atom.
	// Therefor mutate Gly to Ala first
	if(res1.aa() == chemical::aa_gly) {
		//dummy conformation; would only be used if bb atoms missing, e.g. Pro
		conformation::Conformation conformation;
		chemical::ResidueTypeSetCAP restype_set =
			chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );
		res1_ptr = conformation::ResidueCAP(new conformation::Residue(
			restype_set->name_map("ALA"), res1, conformation) );
	}
	if(res2.aa() == chemical::aa_gly) {
		//dummy conformation; would only be used if bb atoms missing, e.g. Pro
		conformation::Conformation conformation;
		chemical::ResidueTypeSetCAP restype_set =
			chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );
		res2_ptr = conformation::ResidueCAP(new conformation::Residue(
			restype_set->name_map("ALA"), res2, conformation) );
	}
	//Make sure they both have CB now
	assert(res1_ptr->type().has("CB"));
	assert(res2_ptr->type().has("CB"));

	Vector const& calpha_1 ( res1_ptr->xyz("CA") );
	Vector const& cbeta_1  ( res1_ptr->xyz("CB") );
	Vector const& n_1      ( res1_ptr->xyz("N")  );
	Vector const& calpha_2 ( res2_ptr->xyz("CA") );
	Vector const& cbeta_2  ( res2_ptr->xyz("CB") );
	Vector const& c_2      ( res2_ptr->xyz("C")  );

	cbcb_distance_sq       = cbeta_1.distance_squared(cbeta_2);
	cacbcb_angle_1      = angle_of( calpha_1, cbeta_1, cbeta_2);
	cacbcb_angle_2      = angle_of( calpha_2, cbeta_2, cbeta_1);
	cacbcb_angle_1 *= radians_to_degrees; // convert
	cacbcb_angle_2 *= radians_to_degrees; // convert
	cacbcbca_dihedral   = dihedral_degrees(calpha_1,cbeta_1,cbeta_2,calpha_2);
	//Use N-Ca-Ca-C instead of N-Ca-Ca-N to follow the fold tree.
	backbone_dihedral   = dihedral_degrees(n_1, calpha_1, calpha_2, c_2);

	centroid_distance_sq = -1;
	if( res1_ptr->type().has("CEN") &&
		res2_ptr->type().has("CEN") )
	{
		Vector const& cen_1( res1_ptr->xyz("CEN"));
		Vector const& cen_2( res2_ptr->xyz("CEN"));
		centroid_distance_sq = cen_1.distance_squared(cen_2);

	//Postcondition validation
		assert(0. <= centroid_distance_sq );
	}
	assert(0. <= cbcb_distance_sq);
	assert(0. <= cacbcb_angle_1); assert( cacbcb_angle_1 <= 180. );
	assert(0. <= cacbcb_angle_2); assert( cacbcb_angle_2 <= 180. );
}

Cb_Distance_FuncCOP CentroidDisulfidePotential::cb_distance_func_ =
	new Cb_Distance_Func();
Cen_Distance_FuncCOP CentroidDisulfidePotential::cen_distance_func_ =
	new Cen_Distance_Func();
CaCbCb_Angle_FuncCOP CentroidDisulfidePotential::cacbcb_angle_func_ =
	new CaCbCb_Angle_Func();
NCaCaC_Dihedral_FuncCOP CentroidDisulfidePotential::ncacac_dihedral_func_ =
	new NCaCaC_Dihedral_Func();
CaCbCbCa_Dihedral_FuncCOP CentroidDisulfidePotential::cacbcbca_dihedral_func_ =
	new CaCbCbCa_Dihedral_Func();

/**
 * @brief Decide whether there is a disulfide bond between two residues.
 *
 * Does not require that the residues be cysteines, so if this is important
 * you should check for CYS first. (The relaxed requirements are useful for
 * design.)
 */
bool CentroidDisulfidePotential::is_disulfide(
		Residue const & res1,
		Residue const & res2) const
{
	//Get Cb distance score
	Energy cbcb_distance_score,
		centroid_distance_score,
		cacbcb_angle_1_score,
		cacbcb_angle_2_score,
		cacbcbca_dihedral_score,
		backbone_dihedral_score;
	Real cbcb_distance_sq,
		centroid_distance_sq,
		cacbcb_angle_1,
		cacbcb_angle_2,
		cacbcbca_dihedral,
		backbone_dihedral,
		score_factor;
	score_disulfide(res1,
		res2,
		cbcb_distance_sq,
		centroid_distance_sq,
		cacbcb_angle_1,
		cacbcb_angle_2,
		cacbcbca_dihedral,
		backbone_dihedral,
		cbcb_distance_score,
		centroid_distance_score,
		cacbcb_angle_1_score,
		cacbcb_angle_2_score,
		cacbcbca_dihedral_score,
		backbone_dihedral_score,
		score_factor
		);


	return cbcb_distance_score <= disulfide_cb_dist_cutoff &&
		cacbcb_angle_1 >= 60. &&
		cacbcb_angle_2 >= 60. ;
}

///@brief the Cysteines with cb dist scores less than this threshold are
/// very likely (99%) to be disulfide bonded.
const Real CentroidDisulfidePotential::disulfide_cb_dist_cutoff(4.392);


///////////////////////
// Scoring Functions //
///////////////////////

/// @brief Helper function for initializing Histograms from the database
/// @note The static functions in FullatomDisulfidePotential are a more elegant
///  way to initialize the Histograms
static HistogramCOP<Real,Real>::Type
histogram_from_db(string file) {
	utility::io::izstream scores_stream;
	basic::database::open( scores_stream, file);
	HistogramCOP<Real,Real>::Type scores = new Histogram<Real,Real>( scores_stream() );
	scores_stream.close();
	return scores;
}

//Cb_Distance_Func
Cb_Distance_Func::Cb_Distance_Func() {}
Cb_Distance_Func::~Cb_Distance_Func() {}
Real Cb_Distance_Func::func( Real const cb_dist_sq_) const {
	using core::scoring::constraints::dgaussian;
	Energy score = base_score_;
	for(Size i = 0; i<3; ++i ) {
		score -= dgaussian(cb_dist_sq_, means_[i], sds_[i], weights_[i] );
	}

	return score;
}
Real Cb_Distance_Func::dfunc( Real const) const {
	return 0.0;
}
const Real Cb_Distance_Func::means_[3] = { 12.445, 15.327, 14.0 };
const Real Cb_Distance_Func::sds_[3]   = { 1.1737973, 2.1955666, 0.3535534 };
const Real Cb_Distance_Func::weights_[3] = {10.8864116, 33.5711622, 0.2658681 };
const Real Cb_Distance_Func::base_score_ = 0.0;

//Cen_Distance_Func
Cen_Distance_Func::Cen_Distance_Func() {}
Cen_Distance_Func::~Cen_Distance_Func() {}
Real Cen_Distance_Func::func( Real const cen_dist_sq) const {
	if( centroid_dist_scores_ == 0)
		centroid_dist_scores_ = histogram_from_db("scoring/score_functions/disulfides/centroid_distance_score");
	Real e(0.0);
	centroid_dist_scores_->interpolate(cen_dist_sq,e);
	return e;
}
Real Cen_Distance_Func::dfunc( Real const ) const {
	return 0.0;
}
HistogramCOP<Real,Real>::Type Cen_Distance_Func::centroid_dist_scores_ = 0;

//CaCbCb_Angle_Func
CaCbCb_Angle_Func::CaCbCb_Angle_Func() {}
CaCbCb_Angle_Func::~CaCbCb_Angle_Func() {}
Real CaCbCb_Angle_Func::func( Real const cacbcb_angle) const {
	if( CaCbCb_angle_scores_ == 0 )
		CaCbCb_angle_scores_ = histogram_from_db("scoring/score_functions/disulfides/centroid_CaCbCb_angle_score");
	Real e(0.0);
	CaCbCb_angle_scores_->interpolate(cacbcb_angle,e);
	return e;
}
Real CaCbCb_Angle_Func::dfunc( Real const ) const {
	return 0.0;
}
HistogramCOP<core::Real,core::Real>::Type CaCbCb_Angle_Func::CaCbCb_angle_scores_ = 0;


//NCaCaC_Dihedral_Func
NCaCaC_Dihedral_Func::NCaCaC_Dihedral_Func() {}
NCaCaC_Dihedral_Func::~NCaCaC_Dihedral_Func() {}
Real NCaCaC_Dihedral_Func::func( Real const backbone_dihedral) const {
	if( backbone_dihedral_scores_ == 0 )
		backbone_dihedral_scores_ = histogram_from_db("scoring/score_functions/disulfides/centroid_backbone_dihedral_score");
	Real e(0.0);
	backbone_dihedral_scores_->interpolate(backbone_dihedral, e);
	return e;
}
Real NCaCaC_Dihedral_Func::dfunc( Real const ) const {
	return 0.0;
}
HistogramCOP<core::Real,core::Real>::Type NCaCaC_Dihedral_Func::backbone_dihedral_scores_ = 0;



//CaCbCbCa_Dihedral_Func
CaCbCbCa_Dihedral_Func::CaCbCbCa_Dihedral_Func() {}
CaCbCbCa_Dihedral_Func::~CaCbCbCa_Dihedral_Func() {}
Real CaCbCbCa_Dihedral_Func::func( Real const cacbcbca_dihedral) const {
	if( CaCbCbCa_dihedral_scores_ == 0 )
		CaCbCbCa_dihedral_scores_ = histogram_from_db("scoring/score_functions/disulfides/centroid_CaCbCbCa_dihedral_score");
	Real e(0.0);
	CaCbCbCa_dihedral_scores_->interpolate(cacbcbca_dihedral,e);
	return e;
}
Real CaCbCbCa_Dihedral_Func::dfunc( Real const ) const {
	return 0.0;
}
HistogramCOP<core::Real,core::Real>::Type CaCbCbCa_Dihedral_Func::CaCbCbCa_dihedral_scores_ = 0;


} // disulfides
} // scoring
} // core
