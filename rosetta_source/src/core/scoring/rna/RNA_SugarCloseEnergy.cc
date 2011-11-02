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
#include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>


// C++


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
	o4star_c1star_bond_length_( 1.414 ),
	o4star_c1star_sd_( 0.01 ),
	o4star_c1star_dist_harm_func_( new constraints::HarmonicFunc( o4star_c1star_bond_length_, scale_rna_torsion_sd_ * o4star_c1star_sd_ )),
	angle_sd_( numeric::conversions::radians( 1.0 ) ),
	o4star_c1star_c2star_bond_angle_( numeric::conversions::radians( 106.39 ) ),
	o4star_c1star_c2star_angle_harm_func_(
		new constraints::HarmonicFunc( o4star_c1star_c2star_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) ),
	o4star_c1star_first_base_bond_angle_( numeric::conversions::radians( 108.2 ) ),
	o4star_c1star_first_base_angle_harm_func_(
		new constraints::HarmonicFunc( o4star_c1star_first_base_bond_angle_, angle_sd_ ) ),
	c4star_o4star_c1star_bond_angle_( numeric::conversions::radians( 110.4 ) ),
	c4star_o4star_c1star_angle_harm_func_(
																				new constraints::HarmonicFunc( c4star_o4star_c1star_bond_angle_, scale_rna_torsion_sd_ * angle_sd_ ) )
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
		pose::Pose const & pose,
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

	if ( !rsd.is_RNA() ) return;

	Size const & i( rsd.seqpos() );

	//String lookups are slow, but needed here to be careful.
	id::AtomID o4star_id( rsd.atom_index( " O4*" ), i );
	id::AtomID c1star_id( rsd.atom_index( " C1*" ), i );
	id::AtomID c2star_id( rsd.atom_index( " C2*" ), i );
	id::AtomID c4star_id( rsd.atom_index( " C4*" ), i );
	id::AtomID first_base_atom_id( first_base_atom_index( rsd ),  i );

	constraints::ConstraintOP dist_cst = new constraints::AtomPairConstraint( o4star_id,
																																						c1star_id,
																																						o4star_c1star_dist_harm_func_,
																																						rna_sugar_close );

	cst_set.add_constraint( dist_cst );
	//	std::cout << "O4*-C1* distance" << rsd.seqpos() << ' ' << o4star_c1star_dist_harm_func_->func( (rsd.xyz(" O4*") - rsd.xyz( " C1*" )).length() ) << std::endl;


	constraints::ConstraintOP angle1 = new constraints::AngleConstraint( o4star_id,
																																			 c1star_id,
																																			 c2star_id,
																																			 o4star_c1star_c2star_angle_harm_func_,
																																			 rna_sugar_close );

	constraints::ConstraintOP angle2 = new constraints::AngleConstraint( c4star_id,
																																			 o4star_id,
																																			 c1star_id,
																																			 c4star_o4star_c1star_angle_harm_func_,
																																			 rna_sugar_close );
	cst_set.add_constraint( angle2 );
	//	std::cout << "C4*-O4*-C1* angle" << rsd.seqpos() <<	 ' '<< c4star_o4star_c1star_angle_harm_func_->func( angle_radians( rsd.xyz( " C4*" ), rsd.xyz( " O4*" ), rsd.xyz( " C1*" ) ) ) << std::endl;


	Size const first_base_index = first_base_atom_index( rsd );
	constraints::ConstraintOP angle3 = new constraints::AngleConstraint( o4star_id,
																																			 c1star_id,
																																			 first_base_atom_id,
																																			 o4star_c1star_first_base_angle_harm_func_,
																																			 rna_sugar_close );

	cst_set.add_constraint( angle3 );
	//	std::cout << "O4*-C1*-FIRSTBASE angle" << rsd.seqpos() <<	 ' ' << o4star_c1star_first_base_angle_harm_func_->func( angle_radians( rsd.xyz( " O4*" ), rsd.xyz( " C1*" ), rsd.xyz( first_base_index ) ) ) << std::endl;

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

