// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJEnergyTable.cc
/// @brief  Molecular mechanics lj energy table class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

// Unit headers
#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMLJScore.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>

// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>

// AUTO-REMOVED #include <basic/prof.hh>

// Utility header
// AUTO-REMOVED #include <utility/keys/Key4Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ headers
// AUTO-REMOVED #include <iostream>
#include <string>
#include <map>
// AUTO-REMOVED #include <math.h>

#include <core/scoring/mm/MMLJLibrary.hh>
#include <utility/vector1.hh>
#include <cmath>


namespace core {
namespace scoring {
namespace mm {

/// @details Auto-generated virtual destructor
MMLJEnergyTable::~MMLJEnergyTable() {}

static basic::Tracer TR("core.scoring.mm.MMLJEnergyTable");

MMLJEnergyTable::MMLJEnergyTable() :
  max_dist_(49.0),
  bin_dist_(0.05),
  bins_per_angstrom_squared_(20),
	linear_switch_point(0.6) // 60% of distance at minimum
{
  using namespace core;
  using namespace chemical;

  // get MMAtomTypeSet associated with the library
  MMAtomTypeSetCAP mm_atom_set( mm_lj_score_.mm_lj_library().mm_atom_set() );
  Size natomtypes( mm_atom_set->n_atomtypes() );

  TR << "Initializing MM LJ Energy Tables with " << natomtypes << " atom types" << std::endl;

  // init both table dimentions for all the tables based on the size of the library
  mm_atom_pair_rep_energy_table_.resize( natomtypes );
  mm_atom_pair_atr_energy_table_.resize( natomtypes );
  mm_atom_pair_rep_deriv_table_.resize( natomtypes );
  mm_atom_pair_atr_deriv_table_.resize( natomtypes );
  mm_atom_pair_rep_three_bond_energy_table_.resize( natomtypes );
  mm_atom_pair_atr_three_bond_energy_table_.resize( natomtypes );
  mm_atom_pair_rep_three_bond_deriv_table_.resize( natomtypes );
  mm_atom_pair_atr_three_bond_deriv_table_.resize( natomtypes );

  for ( Size i = 1; i <= natomtypes; ++i ) {
    mm_atom_pair_rep_energy_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_atr_energy_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_rep_deriv_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_atr_deriv_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_rep_three_bond_energy_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_atr_three_bond_energy_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_rep_three_bond_deriv_table_[i].resize( natomtypes, NULL );
    mm_atom_pair_atr_three_bond_deriv_table_[i].resize( natomtypes, NULL );
  }

	// create MM atom pair energy tables for greater than three bonds
	TR << "Precomputing >=4 bond energy values" << std::endl;
	for ( Size i = 1; i <= natomtypes; ++i ) {
    for ( Size j = 1; j <= natomtypes; ++j ) {
      // if the energy vector at position i,j is not defined check to see if the one at j,i is
      // since they shouold be the same. If the energy vector at j,i is not defined, create it at i,j
      if ( mm_atom_pair_rep_energy_table_[i][j] == NULL ) {
				if( mm_atom_pair_rep_energy_table_[j][i] == NULL ) {
					// create the vector
					mm_atom_pair_rep_energy_table_[i][j] = new EnergyVector;
					mm_atom_pair_atr_energy_table_[i][j] = new EnergyVector;
					mm_atom_pair_rep_deriv_table_[i][j] = new EnergyVector;
					mm_atom_pair_atr_deriv_table_[i][j] = new EnergyVector;

					// get the distance and energy for when the function is a minimum
					Real min_ener_dist( mm_lj_score_.min_dist( i, j, 4 ) );
					Real min_ener( mm_lj_score_.score( i, j, 4, min_ener_dist ) );

					// get values for linear switching at short distances
					Real switch_dist( linear_switch_point * min_ener_dist );
					Real switch_dist_squared( switch_dist * switch_dist );
					Real switch_slope( mm_lj_score_.deriv_score( i, j, 4, switch_dist ) );
					Real switch_ener( mm_lj_score_.score( i, j, 4, switch_dist ) );
					Real switch_intercept( -1 * switch_dist * switch_slope + switch_ener );

					// now the rest, starting from mm_lj_table_bin_dist_ distance
					Real prev_score( switch_intercept );
					for ( Real k = 0; k <= max_dist_; k += bin_dist_ ) {
						Real temp_score(0), temp_deriv(0);
						if( k <= switch_dist_squared ){ // in switch region
							temp_score = switch_slope * sqrt(k) + switch_intercept;
							temp_deriv = switch_slope;
						} else { // not in switch region
							temp_score = mm_lj_score_.score( i, j, 4, sqrt(k) );
							temp_deriv = mm_lj_score_.deriv_score( i, j, 4, sqrt(k) );
						}
						if ( temp_score <= prev_score ) { // repulcive
							mm_atom_pair_rep_energy_table_[i][j]->push_back( temp_score - min_ener );
							mm_atom_pair_atr_energy_table_[i][j]->push_back( min_ener );
							mm_atom_pair_rep_deriv_table_[i][j]->push_back( temp_deriv );
							mm_atom_pair_atr_deriv_table_[i][j]->push_back( 0 );
						} else { // atractive
							mm_atom_pair_rep_energy_table_[i][j]->push_back( 0 );
							mm_atom_pair_atr_energy_table_[i][j]->push_back( temp_score );
							mm_atom_pair_rep_deriv_table_[i][j]->push_back( 0 );
							mm_atom_pair_atr_deriv_table_[i][j]->push_back( temp_deriv );
						}
						prev_score = temp_score;
					}
				} else {
					mm_atom_pair_rep_energy_table_[i][j] = mm_atom_pair_rep_energy_table_[j][i];
					mm_atom_pair_atr_energy_table_[i][j] = mm_atom_pair_atr_energy_table_[j][i];
					mm_atom_pair_rep_deriv_table_[i][j] = mm_atom_pair_rep_deriv_table_[j][i];
					mm_atom_pair_atr_deriv_table_[i][j] = mm_atom_pair_atr_deriv_table_[j][i];
				}
      }
    } // foreach i
  } // foreach j

  // create the atom pair three bond energy tables
	TR << "Precomputing 3 bond energy values" << std::endl;
  for ( Size i = 1; i <= natomtypes; ++i ) {
    for ( Size j = 1; j <= natomtypes; ++j ) {
      // if the energy vector at position i,j is not defined check to see if the one at j,i is
      // since they shouold be the same. If the energy vector at j,i is not defined, create it at i,j
      if ( mm_atom_pair_rep_three_bond_energy_table_[i][j] == NULL ) {
				if( mm_atom_pair_rep_three_bond_energy_table_[j][i] == NULL ) {
					// create the vector
					mm_atom_pair_rep_three_bond_energy_table_[i][j] = new EnergyVector;
					mm_atom_pair_atr_three_bond_energy_table_[i][j] = new EnergyVector;
					mm_atom_pair_rep_three_bond_deriv_table_[i][j] = new EnergyVector;
					mm_atom_pair_atr_three_bond_deriv_table_[i][j] = new EnergyVector;

					// get the distance and energy for when the function is a minimum
					Real min_ener_dist( mm_lj_score_.min_dist( i, j, 3 ) );
					Real min_ener( mm_lj_score_.score( i, j, 3, min_ener_dist ) );

					// get values for linear switching at short distances
					Real switch_dist( linear_switch_point * min_ener_dist );
					Real switch_dist_squared( switch_dist * switch_dist );
					Real switch_slope( mm_lj_score_.deriv_score( i, j, 3, switch_dist ) );
					Real switch_ener( mm_lj_score_.score( i, j, 3, switch_dist ) );
					Real switch_intercept( -1 * switch_dist * switch_slope + switch_ener );

					// now the rest starting from mm_lj_table_bin_dist_ distance
					Real prev_score( switch_intercept );
					for ( Real k = 0; k <= max_dist_; k += bin_dist_ ) {
						Real temp_score(0), temp_deriv(0);
						if( k <= switch_dist_squared ){ // in switch region
							temp_score = switch_slope * sqrt(k) + switch_intercept;
							temp_deriv = switch_slope;
						} else {
							temp_score = mm_lj_score_.score( i, j, 3, sqrt(k) );
							temp_deriv = mm_lj_score_.deriv_score( i, j, 3, sqrt(k) );
						}
						if ( temp_score <= prev_score ) { // repulcive
							mm_atom_pair_rep_three_bond_energy_table_[i][j]->push_back( temp_score - min_ener );
							mm_atom_pair_atr_three_bond_energy_table_[i][j]->push_back( min_ener );
							mm_atom_pair_rep_three_bond_deriv_table_[i][j]->push_back( temp_deriv );
							mm_atom_pair_atr_three_bond_deriv_table_[i][j]->push_back( 0 );
						} else { // atractive
							mm_atom_pair_rep_three_bond_energy_table_[i][j]->push_back( 0 );
							mm_atom_pair_atr_three_bond_energy_table_[i][j]->push_back( temp_score );
							mm_atom_pair_rep_three_bond_deriv_table_[i][j]->push_back( 0 );
							mm_atom_pair_atr_three_bond_deriv_table_[i][j]->push_back( temp_deriv );
						}
						prev_score = temp_score;
					}
				} else {
					mm_atom_pair_rep_three_bond_energy_table_[i][j] = mm_atom_pair_rep_three_bond_energy_table_[j][i];
					mm_atom_pair_atr_three_bond_energy_table_[i][j] = mm_atom_pair_atr_three_bond_energy_table_[j][i];
					mm_atom_pair_rep_three_bond_deriv_table_[i][j] = mm_atom_pair_rep_three_bond_deriv_table_[j][i];
					mm_atom_pair_atr_three_bond_deriv_table_[i][j] = mm_atom_pair_atr_three_bond_deriv_table_[j][i];
				}
      }
    } // foreach i
  } // foreach j
}

void
MMLJEnergyTable::score( Size atom1, Size atom2, Size & path_distance, Real & squared_distance, Real & rep, Real & atr  ) const
{
  // init values
  rep = atr = 0;

  if ( squared_distance <= max_dist_ ) {

    // calc distance bins and frac
    Real bin( squared_distance * bins_per_angstrom_squared_ );
    Size l_bin( static_cast< Size >( bin ) + 1 );
    Size u_bin( l_bin + 1 );
    Real frac( bin - ( l_bin - 1 ) );

		// get correct vectors
		EnergyVector & rep_vec = ( path_distance == 3 ? *mm_atom_pair_rep_three_bond_energy_table_[ atom1][ atom2] : *mm_atom_pair_rep_energy_table_[ atom1][ atom2 ] );
		EnergyVector & atr_vec = ( path_distance == 3 ? *mm_atom_pair_atr_three_bond_energy_table_[ atom1][ atom2] : *mm_atom_pair_atr_energy_table_[ atom1][ atom2 ] );

    // linear interpolate between upper and lower bins
    rep = rep_vec[l_bin] + frac * ( rep_vec[u_bin] - rep_vec[l_bin] );
    atr = atr_vec[l_bin] + frac * ( atr_vec[u_bin] - atr_vec[l_bin] );
  }
}

void
MMLJEnergyTable::deriv_score( Size atom1, Size atom2, Size & path_distance, Real & squared_distance, Real & drep, Real & datr ) const
{
  // init values
  drep = datr = 0;

  if ( squared_distance <= max_dist_ ) {

    // calc distance bins and frac
    Real bin( squared_distance * bins_per_angstrom_squared_ );
    Size l_bin( static_cast< Size >( bin ) + 1 );
    Size u_bin( l_bin + 1 );
    Real frac( bin - ( l_bin - 1 ) );

		// get correct vectors
		EnergyVector & drep_vec = ( path_distance == 3 ? *mm_atom_pair_rep_three_bond_deriv_table_[ atom1][ atom2] : *mm_atom_pair_rep_deriv_table_[ atom1][ atom2 ] );
		EnergyVector & datr_vec = ( path_distance == 3 ? *mm_atom_pair_atr_three_bond_deriv_table_[ atom1][ atom2] : *mm_atom_pair_atr_deriv_table_[ atom1][ atom2 ] );

    // linear interpolate between upper and lower bins
    drep = drep_vec[l_bin] + frac * ( drep_vec[u_bin] - drep_vec[l_bin] );
    datr = datr_vec[l_bin] + frac * ( datr_vec[u_bin] - datr_vec[l_bin] );
  }
}

} // namespace mm
} // namespace scoring
} // namespace core


//  // DEBUG print stuff
//   for ( Size i = 1; i <= natomtypes; ++i ) {
//     for ( Size j = 1; j <= natomtypes; ++j ) {
//       std::cout << "REP_SCORE " << i << ":" << j << "\t" << mm_atom_pair_rep_energy_table_[i][j] << "\t";
//       for ( Size k = 1; k <= mm_atom_pair_rep_energy_table_[i][j]->size(); ++k ) {
// 	std::cout << mm_atom_pair_rep_energy_table_[i][j]->at(k) << " ";
//       }
//       std::cout << std::endl;
//     }
//   }

//   for ( Size i = 1; i <= natomtypes; ++i ) {
//     for ( Size j = 1; j <= natomtypes; ++j ) {
//       std::cout << "ATR_SCORE " << i << ":" << j << "\t" << mm_atom_pair_atr_energy_table_[i][j] << "\t";
//       for ( Size k = 1; k <= mm_atom_pair_atr_energy_table_[i][j]->size(); ++k ) {
// 	std::cout << mm_atom_pair_atr_energy_table_[i][j]->at(k) << " ";
//       }
//       std::cout << std::endl;
//     }
//   }

//   for ( Size i = 1; i <= natomtypes; ++i ) {
//     for ( Size j = 1; j <= natomtypes; ++j ) {
//       std::cout << "REP_DERIV " << i << ":" << j << "\t" << mm_atom_pair_rep_deriv_table_[i][j] << "\t";
//       for ( Size k = 1; k <= mm_atom_pair_rep_deriv_table_[i][j]->size(); ++k ) {
// 	std::cout << mm_atom_pair_rep_deriv_table_[i][j]->at(k) << " ";
//       }
//       std::cout << std::endl;
//     }
//   }

//   for ( Size i = 1; i <= natomtypes; ++i ) {
//     for ( Size j = 1; j <= natomtypes; ++j ) {
//       std::cout << "ATR_DERIV " << i << ":" << j << "\t" << mm_atom_pair_atr_deriv_table_[i][j] << "\t";
//       for ( Size k = 1; k <= mm_atom_pair_atr_deriv_table_[i][j]->size(); ++k ) {
// 	std::cout << mm_atom_pair_atr_deriv_table_[i][j]->at(k) << " ";
//       }
//       std::cout << std::endl;
//     }
//   }

//   EnergyVector* rep_vec;
//   EnergyVector* atr_vec;

//   if( path_distance == 3 ) {
//     rep_vec = mm_atom_pair_rep_three_bond_energy_table_[atom1][atom2];
//     atr_vec = mm_atom_pair_atr_three_bond_energy_table_[atom1][atom2];
//   } else {
//     rep_vec = mm_atom_pair_rep_energy_table_[atom1][atom2];
//     atr_vec = mm_atom_pair_atr_energy_table_[atom1][atom2];
//   }

//   rep = rep_vec->at(l_bin) + frac * ( rep_vec->at(u_bin) - rep_vec->at(l_bin) );
//   atr = atr_vec->at(l_bin) + frac * ( atr_vec->at(u_bin) - atr_vec->at(l_bin) );

//   EnergyVector* drep_vec;
//   EnergyVector* datr_vec;

//   if( path_distance == 3 ) {
//     drep_vec = mm_atom_pair_rep_three_bond_deriv_table_[atom1][atom2];
//     datr_vec = mm_atom_pair_atr_three_bond_deriv_table_[atom1][atom2];
//   } else {
//     drep_vec = mm_atom_pair_rep_deriv_table_[atom1][atom2];
//     datr_vec = mm_atom_pair_atr_deriv_table_[atom1][atom2];
//   }

//   drep = drep_vec->at(l_bin) + frac * ( drep_vec->at(u_bin) - drep_vec->at(l_bin) );
//   datr = datr_vec->at(l_bin) + frac * ( datr_vec->at(u_bin) - datr_vec->at(l_bin) );

