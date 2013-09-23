// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_SugarCloseEnergy.cc
/// @brief  RNA_SugarClose energy method class implementation
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_SugarCloseEnergy.hh>
#include <core/scoring/rna/RNA_SugarCloseEnergyCreator.hh>

// Package Headers
#include <core/chemical/rna/RNA_Util.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/FadeFunc.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>


// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


// C++

using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace rna {

/// @details This must return a fresh instance of the RNA_SugarCloseEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_SugarCloseEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_SugarCloseEnergy;
}

ScoreTypes
RNA_SugarCloseEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_sugar_close );
	return sts;
}


/// ctor
RNA_SugarCloseEnergy::RNA_SugarCloseEnergy() :
	parent( new RNA_SugarCloseEnergyCreator ),
	scale_rna_torsion_tether_( 0.05 ), // THIS IS A SCALING FACTOR FOR ALL CONSTRAINTS.
	scale_rna_torsion_sd_( 1.0 / std::sqrt( scale_rna_torsion_tether_ ) ),
	o4prime_c1prime_bond_length_( 1.414 ),
	o4prime_c1prime_sd_( 0.01 ),
	o4prime_c1prime_dist_harm_func_( new constraints::HarmonicFunc( o4prime_c1prime_bond_length_, scale_rna_torsion_sd_ * o4prime_c1prime_sd_ )),
	angle_sd_( numeric::conversions::radians( 1.0 ) ),
	o4prime_c1prime_c2prime_bond_angle_( numeric::conversions::radians( 106.39 ) ),
	o4prime_c1prime_c2prime_angle_harm_func_(
		new constraints::HarmonicFunc( o4prime_c1prime_c2prime_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) ),
	o4prime_c1prime_first_base_bond_angle_( numeric::conversions::radians( 108.2 ) ),
	o4prime_c1prime_first_base_angle_harm_func_(
		new constraints::HarmonicFunc( o4prime_c1prime_first_base_bond_angle_, angle_sd_ ) ),
	c4prime_o4prime_c1prime_bond_angle_( numeric::conversions::radians( 110.4 ) ),
	c4prime_o4prime_c1prime_angle_harm_func_(
																				new constraints::HarmonicFunc( c4prime_o4prime_c1prime_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) ),
	//phenix_based_sugar_close params
	use_phenix_sugar_close_( basic::options::option[ basic::options::OptionKeys::rna::corrected_geo ]() ),
	o4prime_c1prime_bond_north_(1.412),
	o4prime_c1prime_bond_south_(1.415),
	bond_sd_(0.015),
	o4prime_c1prime_c2prime_angle_north_( numeric::conversions::radians(107.6) ),
	o4prime_c1prime_c2prime_angle_south_( numeric::conversions::radians(105.8) ),
	o4prime_c1prime_n1_9_angle_north_( numeric::conversions::radians(108.5) ),
	o4prime_c1prime_n1_9_angle_south_( numeric::conversions::radians(108.2) ),
	c4prime_o4prime_c1prime_angle_north_( numeric::conversions::radians(109.7) ),
	c4prime_o4prime_c1prime_angle_south_( numeric::conversions::radians(109.9) ),
	angle_sd1_( numeric::conversions::radians(1.0) ),
	angle_sd2_( numeric::conversions::radians(1.5) )
{}

RNA_SugarCloseEnergy::~RNA_SugarCloseEnergy() {}

/// clone
methods::EnergyMethodOP
RNA_SugarCloseEnergy::clone() const
{
	return new RNA_SugarCloseEnergy;
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	setup_sugar_ring_closure_constraints( pose );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::residue_energy(
		conformation::Residue const & rsd,
		EnergyMap & emap 	) const {

	constraints::ConstraintSet residue_rna_sugar_close_constraints;
	add_sugar_ring_closure_constraints( rsd, residue_rna_sugar_close_constraints );
	residue_rna_sugar_close_constraints.eval_intrares_energy( rsd, emap );

}

///////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap & emap 	) const {
	return residue_energy( rsd, emap );
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	//rna_sugar_close_constraints_->eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
	// eval_atom_derivative became deprecated while I was working in another branch -- need to find a replacement?
	rna_sugar_close_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
}


/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::setup_sugar_ring_closure_constraints( pose::Pose & pose ) const{
	rna_sugar_close_constraints_ = constraints::ConstraintSetOP( new constraints::ConstraintSet );

	for (Size i = 1; i <= pose.total_residue(); i++ ){
		add_sugar_ring_closure_constraints( pose.residue( i ), *rna_sugar_close_constraints_ );
	}

}

/////////////////////////////////////////////////////////////////////////////////////////
void
RNA_SugarCloseEnergy::add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const {

	using namespace core::scoring::constraints;

	if ( !rsd.is_RNA() ) return;

	Size const & i( rsd.seqpos() );

	//fast look up!
	Size const o4prime_index=rsd.RNA_type().o4prime_atom_index();
	Size const c1prime_index=rsd.RNA_type().c1prime_atom_index();
	Size const c2prime_index=rsd.RNA_type().c2prime_atom_index();
	Size const c4prime_index=rsd.RNA_type().c4prime_atom_index(); 

	//consistency_check
	if(o4prime_index!=7)  utility_exit_with_message("o4prime_id="+ObjexxFCL::string_of(o4prime_index) + "!=7");
	if(c1prime_index!=10) utility_exit_with_message("c1prime_id="+ObjexxFCL::string_of(c1prime_index) +"!=10");
	if(c2prime_index!=11) utility_exit_with_message("c2prime_id="+ObjexxFCL::string_of(c2prime_index) +"!=11");
	if(c4prime_index!=6)  utility_exit_with_message("c4prime_id="+ObjexxFCL::string_of(c4prime_index) +"!=16");

	id::AtomID const o4prime_id( o4prime_index, i );
	id::AtomID const c1prime_id( c1prime_index, i );
	id::AtomID const c2prime_id( c2prime_index, i );
	id::AtomID const c4prime_id( c4prime_index, i );
	id::AtomID const first_base_atom_id( first_base_atom_index( rsd ),  i );

	constraints::ConstraintOP dist_cst, angle1, angle2, angle3;
	if (use_phenix_sugar_close_) {
		Real const delta = rsd.mainchain_torsion( DELTA );
		RNA_FittedTorsionInfo rna_torsion_fitted_info;
		Real const delta_cutoff = rna_torsion_fitted_info.delta_cutoff();
		if ( delta < delta_cutoff ) { //NORTH
			dist_cst = new AtomPairConstraint( o4prime_id, c1prime_id, 
				new HarmonicFunc( o4prime_c1prime_bond_north_, scale_rna_torsion_sd_ * bond_sd_ ) , rna_sugar_close );
			angle1 = new AngleConstraint( o4prime_id, c1prime_id, c2prime_id,
				new HarmonicFunc( o4prime_c1prime_c2prime_angle_north_, scale_rna_torsion_sd_ * angle_sd1_ ), rna_sugar_close );
			angle2 = new AngleConstraint( c4prime_id, o4prime_id, c1prime_id,
				new HarmonicFunc( c4prime_o4prime_c1prime_angle_north_, scale_rna_torsion_sd_ * angle_sd1_ ), rna_sugar_close );
			angle3 = new AngleConstraint( o4prime_id, c1prime_id, first_base_atom_id,
				new HarmonicFunc( o4prime_c1prime_n1_9_angle_north_, scale_rna_torsion_sd_ * angle_sd2_ ), rna_sugar_close );
		} else { //SOUTH
			dist_cst = new AtomPairConstraint( o4prime_id, c1prime_id, 
				new HarmonicFunc( o4prime_c1prime_bond_south_, scale_rna_torsion_sd_ * bond_sd_ ) , rna_sugar_close );
			angle1 = new AngleConstraint( o4prime_id, c1prime_id, c2prime_id,
				new HarmonicFunc( o4prime_c1prime_c2prime_angle_south_, scale_rna_torsion_sd_ * angle_sd1_ ), rna_sugar_close );
			angle2 = new AngleConstraint( c4prime_id, o4prime_id, c1prime_id,
				new HarmonicFunc( c4prime_o4prime_c1prime_angle_south_, scale_rna_torsion_sd_ * angle_sd1_ ), rna_sugar_close );
			angle3 = new AngleConstraint( o4prime_id, c1prime_id, first_base_atom_id,
				new HarmonicFunc( o4prime_c1prime_n1_9_angle_south_, scale_rna_torsion_sd_ * angle_sd2_ ), rna_sugar_close );
		}
	} else {
		dist_cst = 
			new AtomPairConstraint( o4prime_id, c1prime_id, o4prime_c1prime_dist_harm_func_, rna_sugar_close );
		angle1 = 
			new AngleConstraint( o4prime_id, c1prime_id, c2prime_id, o4prime_c1prime_c2prime_angle_harm_func_, rna_sugar_close );
		angle2 = 
			new AngleConstraint( c4prime_id, o4prime_id, c1prime_id, c4prime_o4prime_c1prime_angle_harm_func_, rna_sugar_close );
		angle3 = 
			new AngleConstraint( o4prime_id, c1prime_id, first_base_atom_id, o4prime_c1prime_first_base_angle_harm_func_, rna_sugar_close );
	}

	cst_set.add_constraint( dist_cst );
	cst_set.add_constraint( angle1 ); //Note to Rhiju (12/25/2011): Previously in Trunk version, angle1 was not added to the cst_set!
	cst_set.add_constraint( angle2 );
	cst_set.add_constraint( angle3 );

}



void RNA_SugarCloseEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const{}

core::Size
RNA_SugarCloseEnergy::version() const
{
	return 1; // A new torsion potential (integration from Das lab branch -- Aug 2011)
}


} // rna
} // scoring
} // core

