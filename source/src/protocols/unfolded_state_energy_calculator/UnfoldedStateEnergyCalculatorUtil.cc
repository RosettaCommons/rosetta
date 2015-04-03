// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/UnfoldedStateEnergyCalculator/UnfoldedStateEnergyCalculatorUtil.hh
/// @brief  Utility functions common to both UnfoldedStateEnergyCalculatorjob distributors
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorUtil.hh>

// Project headers
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ headers
#include <algorithm>
#include <cmath>

static thread_local basic::Tracer TR( "protocols.UnfoldedStateEnergyCalculator.UnfoldedStateEnergyCalculatorUtil" );

namespace protocols {
namespace unfolded_state_energy_calculator {

void
calc_all_averages( utility::vector1< core::scoring::EMapVector > unweighted_energies, core::scoring::EMapVector energy_terms )
{
  using namespace core;
  using namespace core::scoring;
  using namespace utility;

  // calc arithmetic avg, boltzmann weighted avg, or mode of fragment central residue energies
  EMapVector mean;
  EMapVector median;
  EMapVector mode;
  EMapVector boltzmann;

  // for each energy term in the EMapVector
  for ( Size i( 1 ); i <= n_score_types; ++i ) {

    // if the energy term has a non-zero weight
    if ( energy_terms[ ScoreType( i ) ] != 0 ) {

      // create a vector of all the energies for this energy term
      vector1< Real > temp_energies;
      temp_energies.resize( unweighted_energies.size() );

      for ( Size j( 1 ); j <= unweighted_energies.size(); ++j ) {
				temp_energies[ j ] =  unweighted_energies[ j ][ ScoreType( i ) ];
      }

      // sort array
      std::sort( temp_energies.begin(), temp_energies.end() );

      // add values to correct EMapVectors
      mean[ ScoreType( i ) ] = calc_vector_mean( temp_energies );
      median[ ScoreType( i ) ] = calc_vector_median( temp_energies );
      mode[ ScoreType( i ) ] = calc_vector_mode( temp_energies );
      boltzmann[ ScoreType( i ) ] = calc_vector_boltzmann( temp_energies );
    }
  }

  // output all the data we have been collecting. despite function name these energies are unweighted
  // since the non-zero weights were set to 1 in set_energy_terms()
	TR << "NUM FRAGMENTS: " << unweighted_energies.size() << std::endl;
	for ( Size i(1); i <= unweighted_energies.size(); ++i ) {
		TR << "FRAG " << i << ":\t" << unweighted_energies[ i ].weighted_string_of( energy_terms ) << std::endl;
	}

	// ouput averaged vectors
  TR << "MEAN UNFODLED ENERGIES:     " << mean.weighted_string_of( energy_terms ) << std::endl;
  TR << "MEDIAN UNFODLED ENERGIES:   " << median.weighted_string_of( energy_terms ) << std::endl;
  TR << "MODE UNFODLED ENERGIES:     " << mode.weighted_string_of( energy_terms ) << std::endl;
  TR << "BOLZMANN UNFOLDED ENERGIES: " << boltzmann.weighted_string_of( energy_terms ) << std::endl;
}


/// @brief assumes sorted array
core::Real
calc_vector_mean( utility::vector1< core::Real> const & data )
{
  using namespace core;

  // sum elements of array
  Real sum( 0 );
  for ( Size i( 1 ); i <= data.size(); ++i ) {
    sum += data[ i ];
  }

  // return average
  return sum/data.size();
}

/// @brief assumes sorted array
core::Real
calc_vector_median(  utility::vector1< core::Real> const & /*data*/ )
{
  using namespace core;

  return 0.0;
}

/// @brief assumes sorted array
core::Real
calc_vector_mode( utility::vector1< core::Real> const & /*data*/ )
{
  using namespace core;

  // calc mode

  // if no mode try reducing signifigance

  return 0.0;
}

/// @brief assumes sorted array
core::Real
calc_vector_boltzmann( utility::vector1< core::Real> const & data )
{
  using namespace core;
  using namespace utility;

  // hard coded temp (K) and boltzmann constant (kcals mol^-1 K^-1)
  Real T( 300.00 ), R( 0.0019872 );

  // calc energy differences
  Real low_energy( data[ 1 ] );
  vector1< core::Real> energy_diffs;
  energy_diffs.resize( data.size() );
  for ( Size i( 1 ); i <= data.size(); ++i ) {
    energy_diffs[ i ] =  data[ i ] - low_energy;
  }

  // calc boltzmann factors
  vector1< core::Real> bolzmann_factors;
  bolzmann_factors.resize( energy_diffs.size() );
  for ( Size i( 1 ); i <= energy_diffs.size(); ++i ) {
    bolzmann_factors[ i ] = exp( ( -1 * energy_diffs[ i ] ) / ( R * T ) );
  }

  // calc average boltzmann factor
  Real sum( 0 ), avg( 0 );
  for ( Size i( 1 ); i <= bolzmann_factors.size(); ++i ) {
    sum += bolzmann_factors[ i ];
  }
  avg = sum/bolzmann_factors.size();

  // return boltzmann weighted average energy
  return (-1 * R * T * log( avg ) ) + low_energy;
}

} // UnfoldedStateEnergyCalculator
} // protocols
