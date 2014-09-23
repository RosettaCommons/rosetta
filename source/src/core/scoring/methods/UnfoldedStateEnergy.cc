// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/UnfoldedStateEnergy.cc
/// @brief  Unfolded state energy method implementation; energies based on eneriges of residues in fragments
/// @author Ron Jacak (ronj@email.unc.edu)

// Unit headers
#include <core/scoring/methods/UnfoldedStateEnergy.hh>
#include <core/scoring/methods/UnfoldedStateEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/UnfoldedStatePotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "core.scoring.methods.UnfoldedStateEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the UnfoldedStateEnergy class, never an instance already in use
methods::EnergyMethodOP
UnfoldedStateEnergyCreator::create_energy_method( methods::EnergyMethodOptions const & options ) const {

	if ( options.has_method_weights( unfolded ) ) {

		utility::vector1< Real > const & v = options.method_weights( unfolded );
		assert( v.size() == scoring::n_score_types );

		// convert the vector of Reals into an EnergyMap, because that's what the constructor for USE takes.
		// assumes that the vector of Reals coming in contains the weights for each of the score types in the
		// scoring namespace enumeration, and more importantly, in the same order.
		EnergyMap e;
		for ( Size ii=1; ii < scoring::n_score_types; ++ii ) {
			e[ (ScoreType) ii ] = v[ii];
		}

		return new UnfoldedStateEnergy( options.unfolded_energies_type(), e );
	}
	return new UnfoldedStateEnergy( options.unfolded_energies_type() );
}

ScoreTypes
UnfoldedStateEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( unfolded );
	return sts;
}


///
/// UnfoldedStateEnergy class methods
///

UnfoldedStateEnergy::UnfoldedStateEnergy( std::string const & type ) :
	parent( methods::EnergyMethodCreatorOP( new UnfoldedStateEnergyCreator ) ),
	type_( type ),
	unf_state_potential_( ScoringManager::get_instance()->get_UnfoldedStatePotential( type ) ),
	score_type_weights_( unf_state_potential_.get_unfoled_potential_file_weights() )
{
	TR.Debug << "instantiating class with weights: " << score_type_weights_.show_nonzero() << std::endl;
}

UnfoldedStateEnergy::UnfoldedStateEnergy( std::string const & type, const EnergyMap & emap_in ):
	parent( methods::EnergyMethodCreatorOP( new UnfoldedStateEnergyCreator ) ),
	type_( type ),
	unf_state_potential_( ScoringManager::get_instance()->get_UnfoldedStatePotential( type ) ),
	score_type_weights_( emap_in )
{
	TR.Debug << "instantiating class with weights: " << score_type_weights_.show_nonzero()  << std::endl;
}

UnfoldedStateEnergy::~UnfoldedStateEnergy() {}

EnergyMethodOP
UnfoldedStateEnergy::clone() const {
	return new UnfoldedStateEnergy( type_, score_type_weights_ );
}

void
UnfoldedStateEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{

	// the value this function returns depends on how this class was constructed. if a set of weights was
	// passed in when the class was constructed, then this function will return a non-zero value.
	// if no weights were passed in, the member variable energy map will be just zeros and the function
	// will return 0.0. that's what we want during the first phase of optE. once a weight set has been
	// established, this class will get recreated with that set of weights and the function will return
	// meaningful energies.

	EnergyMap unweighted_unfolded_energies;
	unf_state_potential_.raw_unfolded_state_energymap( rsd.type().name3(), unweighted_unfolded_energies );

	// don't forget to include the weight for the entire term. no, wait, that weight should get applied later
	// in the scoring process. this method should just return the unweighted unfolded state energy.
	emap[ unfolded ] += unweighted_unfolded_energies.dot( score_type_weights_ );

	return;
}

void
UnfoldedStateEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}
core::Size
UnfoldedStateEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

