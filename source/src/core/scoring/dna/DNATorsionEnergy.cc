// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DNATorsionEnergy.cc
/// @brief  DNATorsion energy method class implementation
/// @author Phil Bradley
/// @author Mike Tyka (mtyka@u.washington.edu)
/// @author Rhiju Das
/// @author Jim Havranek

// Unit Headers
#include <core/scoring/dna/DNATorsionEnergy.hh>
#include <core/scoring/dna/DNATorsionEnergyCreator.hh>

// Package Headers
//#include <core/scoring/dna/DNATorsionPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


//#include <core/pack/task/PackerTask.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>

// Utility headers
#include <numeric/conversions.hh>

// C++


namespace core {
namespace scoring {
namespace dna {

methods::EnergyMethodOP
DNATorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new DNATorsionEnergy;
}

ScoreTypes
DNATorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_bb_torsion );
	sts.push_back( dna_sugar_close );
	sts.push_back( dna_base_distance );
	return sts;
}


/// ctor
DNATorsionEnergy::DNATorsionEnergy():
	parent( methods::EnergyMethodCreatorOP( new DNATorsionEnergyCreator ) ),
	dna_torsion_potential_( ScoringManager::get_instance()->get_DNATorsionPotential() ),
	constraints_ready_( false ),
	verbose_( false )
{}

/// clone
methods::EnergyMethodOP
DNATorsionEnergy::clone() const
{
	return new DNATorsionEnergy;
}


///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	//	if ( constraints_ready_) return;

	dna_torsion_potential_.setup_constraints( pose, dna_torsion_constraints_, dna_sugar_close_constraints_, dna_base_distance_constraints_ );

	constraints_ready_ = true;
}

///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	//	if ( constraints_ready_) return;
	dna_torsion_potential_.setup_constraints( pose, dna_torsion_constraints_, dna_sugar_close_constraints_, dna_base_distance_constraints_ );
	constraints_ready_ = true;
}

///////////////////////////////////////////////////////////////////////////////
/*
void
DNATorsionEnergy::setup_for_packing( pose::Pose & pose, pack::task::PackerTask const & task ) const{
	//	if ( constraints_ready_) return;
	dna_torsion_potential_.setup_constraints( pose, dna_torsion_constraints_, dna_sugar_close_constraints_, dna_base_distance_constraints_ );

	//If designing, can't implement sugar closure constraints, because one of them depends on the first base atom ...
	// whose name and location differ between bases! Luckily, sugar atoms and first base atom should not vary during design.
	// So clear these constraints...
	if ( task.design_any() ) dna_sugar_close_constraints_ = constraints::ConstraintSetOP( new constraints::ConstraintSet );

	constraints_ready_ = true;
}
*/
///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap ) const {

	if ( dna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	EnergyMap emap_temp;
	dna_torsion_constraints_->residue_pair_energy( rsd1, rsd2, pose, sfxn, emap_temp );
	emap[ dna_bb_torsion] += emap_temp[ dna_bb_torsion];

	EnergyMap emap_2nd;
	dna_base_distance_constraints_->residue_pair_energy( rsd1, rsd2, pose, sfxn, emap_2nd );
	emap[ dna_base_distance ] += emap_2nd[ dna_base_distance ];
//	std::cout << "respair_energy" << rsd1.seqpos() << " " << rsd2.seqpos() <<  " " <<  emap_2nd[ dna_base_distance ] << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap 	) const {

	if ( dna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	////////////////////////////////////////////////////////////////////////////////////
	dna_torsion_constraints_->eval_intrares_energy( rsd, pose, sfxn, emap );
	if ( verbose_ ) std::cout << "accumulated dna_torsion " << rsd.seqpos() << " " <<  emap[ dna_bb_torsion ] << std::endl;

	////////////////////////////////////////////////////////////////////////////////////
	// This seems a little wasteful ... there is a "global" dna_sugar_close_constraints_
	// constraint set. But if the residue changes -- during design for example -- then we need to be more careful,
	// as atom names and indices change around.
	constraints::ConstraintSet residue_dna_sugar_close_constraints;
	dna_torsion_potential_.add_sugar_ring_closure_constraints( rsd, residue_dna_sugar_close_constraints );
	residue_dna_sugar_close_constraints.eval_intrares_energy( rsd, pose, sfxn, emap );
	if ( verbose_ ) std::cout << "accumulated dna_sugar_close " <<   rsd.seqpos() << " " << emap[ dna_sugar_close ] << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const & sfxn,
		EnergyMap & totals
) const
{

	if ( dna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	if ( verbose_ ) std::cout << " final_energy " << totals[ dna_bb_torsion ] << std::endl;
	dna_torsion_constraints_->eval_non_residue_pair_energy( pose, sfxn, totals );
	if ( verbose_ ) std::cout << " final_energy " << totals[ dna_bb_torsion ] << std::endl;

	// Sugar constraints are always intra-residue
	//dna_sugar_close_constraints_->eval_non_residue_pair_energy( pose, sfxn, totals );

	constraints_ready_ = false;

}

///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( dna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;
	dna_torsion_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
	dna_sugar_close_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
	dna_base_distance_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
}
///////////////////////////////////////////////////////////////////////////////
/// side chain torsion tethers...
Real
DNATorsionEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap const &
) const
{
	//	if ( dna_torsion_constraints_ == 0 ) return 0.0;
	//	if ( dna_side_chain_torsion_tethers_ == 0 ) return 0.0;
	if ( !constraints_ready_ ) return 0.0;

	// DOF_constraint machinery ... should remove this call, and perhaps all
	// calls to dof_constraint!!
	//	Real result ( dna_torsion_constraints_->eval_dof_derivative( id, tor, pose, scorefxn, weights ) );

	Real result( 0.0 );

	return result;
}



///////////////////////////////////////////////////////////////////////////////
void
DNATorsionEnergy::finalize_after_derivatives(
		pose::Pose &,
		ScoreFunction const &
) const
{

	if ( dna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;
	constraints_ready_ = false;

}


void DNATorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required */) const {}

core::Size
DNATorsionEnergy::version() const
{
	return 1;
}


/// @brief DNA_PairwiseLowResolutionEnergy distance cutoff
Distance
DNATorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}


} // dna
} // scoring
} // core

