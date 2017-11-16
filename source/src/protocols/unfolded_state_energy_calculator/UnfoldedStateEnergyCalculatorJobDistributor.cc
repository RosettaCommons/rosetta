// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorJobDistributor.cc
/// @brief  Job distributor for UnfoldedStateEnergyCalculator
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers

#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorJobDistributor.hh>
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorUtil.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>

#include <basic/options/keys/OptionKeys.hh>


static basic::Tracer TR( "protocols.UnfoldedStateEnergyCalculator.UnfoldedStateEnergyCalculatorJobDistributor" );

namespace protocols {
namespace unfolded_state_energy_calculator {

using namespace basic::options;
using namespace basic::options::OptionKeys;

/// @brief ctor
UnfoldedStateEnergyCalculatorJobDistributor::UnfoldedStateEnergyCalculatorJobDistributor() :
	FileSystemJobDistributor()
{}

/// @brief dtor (don't put anything in here)
UnfoldedStateEnergyCalculatorJobDistributor::~UnfoldedStateEnergyCalculatorJobDistributor()
= default;

/// @brief
void
UnfoldedStateEnergyCalculatorJobDistributor::go( protocols::moves::MoverOP mover )
{
	using namespace core;
	using namespace core::scoring;
	using namespace utility;

	// call main jd function
	go_main( mover );

	// calc average unweighted energies for all amino acids in the map
	for ( auto & i : unweighted_energies_map_ ) {
		TR << "Calculating averages for " << i.first << std::endl;
		calc_all_averages( i.second, energy_terms_ );
	}
}

void
UnfoldedStateEnergyCalculatorJobDistributor::add_unfolded_energy_data( std::string tlc, core::scoring::EMapVector const & scores )
{
	unweighted_energies_map_[tlc].push_back( scores );
}

/// @details Set the the internal EMapVector that contains the terms in the energy function used to score the
/// fragments. Also if a term has a non-zero weight, set the weight to 1. This allows us to use the EMapVector
/// output weighted functions.
void
UnfoldedStateEnergyCalculatorJobDistributor::set_energy_terms(core::scoring::EMapVector const & weights )
{
	using namespace core;
	using namespace core::scoring;

	// get weights
	energy_terms_ = weights;

	// for each energy term in the EMapVector
	for ( Size i( 1 ); i <= n_score_types; ++i ) {

		// if the energy term has a non-zero weight set it to one
		if ( energy_terms_[ ScoreType( i ) ] > 0 ) {
			energy_terms_[ ScoreType( i ) ] = 1;
		} else if ( energy_terms_[ ScoreType( i ) ] < 0 ) {
			energy_terms_[ ScoreType( i ) ] = -1;
		}
	}
}

} // UnfoldedStateEnergyCalculator
} // protocols
