// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_TorsionEnergy.cc
/// @brief  RNA_Torsion energy method class implementation
/// @author Phil Bradley
/// @author Mike Tyka (mtyka@u.washington.edu)
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/rna/RNA_TorsionEnergy.hh>
#include <core/scoring/rna/RNA_TorsionEnergyCreator.hh>

// Package Headers
//#include <core/scoring/rna/RNA_TorsionPotential.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>



// Project headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
//#include <core/pack/task/PackerTask.hh>

// Utility headers
#include <numeric/conversions.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyMap.hh>


// C++


namespace core {
namespace scoring {
namespace rna {


/// @details This must return a fresh instance of the RNA_TorsionEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_TorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RNA_TorsionEnergy;
}

ScoreTypes
RNA_TorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_torsion );
	return sts;
}


/// ctor
RNA_TorsionEnergy::RNA_TorsionEnergy() :
	parent( new RNA_TorsionEnergyCreator ),
	rna_torsion_potential_( ScoringManager::get_instance()->get_RNA_TorsionPotential() ),
	constraints_ready_( false ),
	verbose_( false )
{}

/// clone
methods::EnergyMethodOP
RNA_TorsionEnergy::clone() const
{
	return new RNA_TorsionEnergy;
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	//	if ( constraints_ready_) return;

	rna_torsion_potential_.setup_constraints( pose, rna_torsion_constraints_, rna_sugar_close_constraints_, rna_side_chain_torsion_tethers_ );

	constraints_ready_ = true;
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	//	if ( constraints_ready_) return;
	rna_torsion_potential_.setup_constraints( pose, rna_torsion_constraints_, rna_sugar_close_constraints_, rna_side_chain_torsion_tethers_ );
	constraints_ready_ = true;
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & designing_residues
) const
{
	//	if ( constraints_ready_) return;
	rna_torsion_potential_.setup_constraints( pose, rna_torsion_constraints_, rna_sugar_close_constraints_, rna_side_chain_torsion_tethers_ );

	//If designing, can't implement sugar closure constraints, because one of them depends on the first base atom ...
	// whose name and location differ between bases! Luckily, sugar atoms and first base atom should not vary during design.
	// So clear these constraints...
	bool designing_any = false;
	for ( Size ii = 1; ii <= designing_residues.size(); ++ii ) {
		if ( designing_residues[ ii ] ) {
			designing_any = true;
			break;
		}
	}
	if ( designing_any ) rna_sugar_close_constraints_ = constraints::ConstraintSetOP( new constraints::ConstraintSet );

	constraints_ready_ = true;
}

///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{

	if ( rna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	EnergyMap emap_temp;
	rna_torsion_constraints_->residue_pair_energy( rsd1, rsd2, pose, sfxn, emap_temp );
	emap[ rna_torsion] += emap_temp[ rna_torsion];

	//if ( verbose_ ) std::cout << "respair_energy" << rsd1.seqpos() << " " << rsd2.seqpos() <<  " " <<  emap_temp[ rna_torsion ] << std::endl;

	// Sugar constraints are always intra-residue
	//rna_sugar_close_constraints_->residue_pair_energy( rsd1, rsd2, pose, sfxn, emap_temp );
	//	emap[ rna_sugar_close] += emap_temp[ rna_sugar_close];

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap 	) const
{

	if ( rna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	////////////////////////////////////////////////////////////////////////////////////
	rna_torsion_constraints_->eval_intrares_energy( rsd, pose, sfxn, emap );
	if ( verbose_ ) std::cout << "accumulated rna_torsion " << rsd.seqpos() << " " <<  emap[ rna_torsion ] << std::endl;

	////////////////////////////////////////////////////////////////////////////////////
	// This seems a little wasteful ... there is a "global" rna_sugar_close_constraints_
	// constraint set. But if the residue changes -- during design for example -- then we need to be more careful,
	// as atom names and indices change around.
	constraints::ConstraintSet residue_rna_sugar_close_constraints;
	rna_torsion_potential_.add_sugar_ring_closure_constraints( rsd, residue_rna_sugar_close_constraints );
	residue_rna_sugar_close_constraints.eval_intrares_energy( rsd, pose, sfxn, emap );
	if ( verbose_ ) std::cout << "accumulated rna_sugar_close " <<   rsd.seqpos() << " " << emap[ rna_sugar_close ] << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////
	// MOVE THIS INSIDE RNA_TorsionPotential, or at least to another function in here.
	// Tethers are assumed to be "intrares", e.g., chi angles.
	Size const i( rsd.seqpos() );

	for (Size n = 1; n <= rna_side_chain_torsion_tethers_[ i ].size(); n++ ) {

		id::TorsionID       const & torsion_id( rna_side_chain_torsion_tethers_[ i ][ n ].first );
		assert( torsion_id.type() == id::CHI );
		Size const chino = torsion_id.torsion();

		utility::vector1< constraints::FuncOP > const & funcs ( rna_side_chain_torsion_tethers_[ i ][ n ].second );
		Real best_value( 0.0 );
		for ( Size m = 1; m <= funcs.size(); m++ ) {
			constraints::FuncOP const & func(  funcs[ m ] );
			Real const value =  func->func( numeric::conversions::radians( rsd.chi( chino ) ) );
			if ( value < best_value || m == 1) best_value = value;
		}
		if (verbose_) std::cout << "sidechain torsion " << i << " " << n << " " << best_value << std::endl;
		emap[ rna_torsion ] += best_value;
	}

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	EnergyMap & totals
) const
{

	if ( rna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;

	if ( verbose_ ) std::cout << " final_energy " << totals[ rna_torsion ] << std::endl;
	rna_torsion_constraints_->eval_non_residue_pair_energy( pose, sfxn, totals );
	if ( verbose_ ) std::cout << " final_energy " << totals[ rna_torsion ] << std::endl;

	// Sugar constraints are always intra-residue
	//rna_sugar_close_constraints_->eval_non_residue_pair_energy( pose, sfxn, totals );

	constraints_ready_ = false;

}

///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( rna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;
	rna_torsion_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
	rna_sugar_close_constraints_->deprecated_eval_atom_derivative( id, pose, sfxn, weights, F1, F2 );
}
///////////////////////////////////////////////////////////////////////////////
/// side chain torsion tethers...
Real
RNA_TorsionEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const & torsion_id_in,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap const & weights
) const
{
	//	if ( rna_torsion_constraints_ == 0 ) return 0.0;
	//	if ( rna_side_chain_torsion_tethers_ == 0 ) return 0.0;
	if ( !constraints_ready_ ) return 0.0;

	// DOF_constraint machinery ... should remove this call, and perhaps all
	// calls to dof_constraint!!
	//	Real result ( rna_torsion_constraints_->eval_dof_derivative( id, tor, pose, scorefxn, weights ) );

	Real result( 0.0 );
	Size const i ( torsion_id_in.rsd() );

	for (Size n = 1; n <= rna_side_chain_torsion_tethers_[ i ].size(); n++ ) {

		id::TorsionID       const & torsion_id( rna_side_chain_torsion_tethers_[ i ][ n ].first );

		if ( torsion_id == torsion_id_in ) {
			assert( torsion_id.type() == id::CHI );
			/*Size const & chino = torsion_id.torsion();*/
			Real const & chi = numeric::conversions::radians( pose.torsion( torsion_id ) );

			utility::vector1< constraints::FuncOP > const & funcs ( rna_side_chain_torsion_tethers_[ i ][ n ].second );
			Real best_value( 0.0 );
			for ( Size m = 1; m <= funcs.size(); m++ ) {
				constraints::FuncOP const & func(  funcs[ m ] );
				Real const value =  func->func( chi );
				if ( value < best_value || m == 1) {
					best_value = value;
					result = weights[ rna_torsion ] * ( func->dfunc( chi ) );
				}
			}
		}

	}

	return result;
}



///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::finalize_after_derivatives(
		pose::Pose &,
		ScoreFunction const &
) const
{

	if ( rna_torsion_constraints_ == 0 ) return;
	if ( !constraints_ready_ ) return;
	constraints_ready_ = false;

}


void RNA_TorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required */) const {}

/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNA_TorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}
core::Size
RNA_TorsionEnergy::version() const
{
	return 1; // Initial versioning
}



} // rna
} // scoring
} // core

