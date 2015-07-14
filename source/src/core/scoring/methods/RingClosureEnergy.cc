// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RingClosureEnergy.cc
/// @brief  Noncanonical ring closure energy method class implementation
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker laboratory

// Unit Headers
#include <core/scoring/methods/RingClosureEnergy.hh>
#include <core/scoring/methods/RingClosureEnergyCreator.hh>

// Package Headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <basic/Tracer.hh>

//C++ header
#include <stdio.h>

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the RingClosureEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RingClosureEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RingClosureEnergy );
}

ScoreTypes
RingClosureEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ring_close );
	return sts;
}


/// @brief Constructor.
///
RingClosureEnergy::RingClosureEnergy() :
	parent( methods::EnergyMethodCreatorOP( new RingClosureEnergyCreator ) ),
	std_dev_sq_(basic::options::option[ basic::options::OptionKeys::score::ring_close_shadow_constraint ])
{
	std_dev_sq_ = std_dev_sq_*std_dev_sq_; //Since we'll only ever use the square of this value, let's store the square of the value, calculated once, rather than recalculating it a zillion times.
}

/// @brief Copy constructor.
///
RingClosureEnergy::RingClosureEnergy( RingClosureEnergy const &src ):
	parent( methods::EnergyMethodCreatorOP( new RingClosureEnergyCreator ) ),
	std_dev_sq_( src.std_dev_sq_ )
{}

/// @brief Clone -- creates a copy and returns an owning pointer to the copy.
///
EnergyMethodOP
RingClosureEnergy::clone() const
{
	return EnergyMethodOP( new RingClosureEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
RingClosureEnergy::residue_energy(
	conformation::Residue const &rsd,
	pose::Pose const & ,//pose,
	EnergyMap &emap
) const
{
	if( rsd.is_virtual_residue() ) return; //Skip virtual residues.
	if( !rsd.has_shadow_atoms() ) return; //Skip residues without shadow atoms.

	for(core::Size ia=1, iamax=rsd.natoms(); ia<=iamax; ++ia) { //Loop through all atoms in the residue.
		core::Size const ia2( rsd.type().atom_being_shadowed( ia ) );
		if(ia2==0) continue; //If this atom doesn't shadow anything, continue.
		
		Distance const distsq ( rsd.xyz( ia ).distance_squared( rsd.xyz( ia2 ) ) ); //Measure the square of the distance between the atom that shadows and the atom being shadowed.		
		emap[ ring_close ] += distsq / ( std_dev_sq_ ); //Note that std_dev_sq_ is actually the SQUARE of the standard deviation.  This is a harmonic potential.	
	}	
	
	return;
}

/// @brief Evaluate the derivatives for all atoms in this residue.
///
void
RingClosureEnergy::eval_residue_derivatives(
	conformation::Residue const &rsd,
	ResSingleMinimizationData const & /*min_data*/,
	pose::Pose const & /*pose*/,
	EnergyMap const &weights,
	utility::vector1< DerivVectorPair > &atom_derivs
) const {

	if( rsd.is_virtual_residue() ) return; //Skip virtual residues.
	if( !rsd.has_shadow_atoms() ) return; //Skip residues without shadow atoms.

	for(core::Size ia=1, iamax=rsd.natoms(); ia<=iamax; ++ia) { //Loop through all atoms in the residue.
		core::Size const ia2( rsd.type().atom_being_shadowed( ia ) );
		if(ia2==0) continue; //If this atom doesn't shadow anything, continue.
		
		//Positions of atom 1 and atom 2:
		Vector const atom1_pos( rsd.xyz( ia  ) );
		Vector const atom2_pos( rsd.xyz( ia2 ) );

		//Storage for f1 and f2 components of the atomic derivatives:
		Vector f1( 0.0 ), f2( 0.0 );
	
		//Storage for the interatomic distance:
		Distance dist( 0.0 );
		
		//Calculate dist, f1, and f2 given atom1_pos and atom2_pos.  Note that f1 and f2 must subsequently be multiplied by the derivative of
		//f(dist) with respect to distance:
		numeric::deriv::distance_f1_f2_deriv( atom1_pos, atom2_pos, dist, f1, f2 );		
		
		//Calculate df(dist)/d_dist.
		//Since f = dist^2/stdev^2, df/d_dist = 2*dist/stdev^2, which must be premultiplied by the ring_close weight:
		core::Real const deriv( weights[ ring_close ] * 2 * dist / std_dev_sq_ );
		
		//Multiply the f1 and f2 vectors by the derivative:
		f1 *= deriv;
		f2 *= deriv;
		
		//Add the f1 and f2 vectors to the accumulators for the respective atoms.
		atom_derivs[ ia  ].f1() += f1;
		atom_derivs[ ia  ].f2() += f2;
		atom_derivs[ ia2 ].f1() -= f1;
		atom_derivs[ ia2 ].f2() -= f2;
		
	} //End loop through atoms.	
	
	return;
}

/// @brief RingClosure Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
RingClosureEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
RingClosureEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

