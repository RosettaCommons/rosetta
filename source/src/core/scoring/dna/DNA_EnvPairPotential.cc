// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/dna/DNA_EnvPairPotential.cc
/// @brief  dna scoring
/// @author Phil Bradley


// Unit Headers
#include <core/scoring/dna/DNA_EnvPairPotential.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

// Package headers

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>

namespace core {
namespace scoring {
namespace dna {

chemical::AA
dna_aa_from_oneletter_code( char const c )
{
	using namespace chemical;
	switch ( c ) {
	case 'a' : return na_ade;
	case 'c' : return na_cyt;
	case 'g' : return na_gua;
	case 't' : return na_thy;
	default :
		utility_exit_with_message("unrecognized oneletter code for dna "+std::string(1,c));
	}
	return aa_unk;
}

Real const DNA_EnvPairPotential::nbr_dis2_threshold_( 12.0 * 12.0 );

Real const DNA_EnvPairPotential::min_env_enrichment_( 0.3333 );
Real const DNA_EnvPairPotential::max_env_enrichment_( 3.0 );
Real const DNA_EnvPairPotential::min_pair_enrichment_( 0.5 );
Real const DNA_EnvPairPotential::max_pair_enrichment_( 2.0 );


/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
DNA_EnvPairPotential::DNA_EnvPairPotential()
{
	// read the database file
	load_score_tables();
}


/// This is really slow. Could be made faster in many ways, but mostly by adding an explicit centroid atom to centroiddna
Vector
DNA_EnvPairPotential::centroid_xyz( conformation::Residue const & rsd ) const
{
	using namespace chemical;
	if ( rsd.is_DNA() ) {
		switch( rsd.aa() ) {
		case na_ade :
			return Real( 0.5 ) * ( rsd.xyz( "N7"  ) + rsd.xyz( "N6" ) );
		case na_gua :
			return Real( 0.5 ) * ( rsd.xyz( "N7"  ) + rsd.xyz( "O6" ) );
		case na_thy :
			return Real( 0.5 ) * ( rsd.xyz( "C7" ) + rsd.xyz( "O4" ) ); // C7 was C5M
		case na_cyt :
			return Real( 0.5 ) * ( rsd.xyz( "C5"  ) + rsd.xyz( "N4" ) );
		default :
			utility_exit_with_message("bad dna aa");
		}
	}
	return rsd.nbr_atom_xyz();
}

Real
DNA_EnvPairPotential::residue_pair_score( AA const & na, AA const & aa, Real const dis2 ) const
{
	using namespace chemical;
	debug_assert( aa <= num_canonical_aas && na >= first_DNA_aa && na <= last_DNA_aa );

	return pair_stats_[ get_pair_disbin( dis2 ) ][ na - first_DNA_aa + 1 ][ aa ];
}

Real
DNA_EnvPairPotential::residue_env_score( AA const & aa, Size const nbr_count ) const
{
	using namespace chemical;
	debug_assert( aa <= num_canonical_aas );

	return env_stats_[ get_env_nbr_bin( nbr_count ) ][ aa ];
}

void
DNA_EnvPairPotential::load_score_tables()
{
	using namespace chemical;

	utility::vector1< Real > twenty_zeroes( 20, 0.0 );
	utility::vector1< utility::vector1< Real > > eighty_zeroes( 4, twenty_zeroes );

	utility::io::izstream data;
	basic::database::open( data, "scoring/dna/dna_env_pair.dat" );


	std::string line;

	while ( getline( data,line ) ) {
		std::istringstream l( line );

		std::string tag;
		char aa('?');
		Real enrich, tmp;
		Size nbr_bin(0);

		l >> tag;
		if ( l.fail() ) continue;
		if ( tag.empty() || tag[0] == '#' ) continue;

		if ( tag == "MAX_PAIR_DIS" ) {
			l >> tmp;
			debug_assert( std::abs( tmp * tmp - nbr_dis2_threshold_ ) < 1e-3 );

		} else if ( tag == "N_PAIR_DISBINS" ) {
			l >> n_pair_disbins_;
			pair_stats_.resize( n_pair_disbins_, eighty_zeroes );
			pair_bin_lower_bounds_.resize( n_pair_disbins_ );

		} else if ( tag == "BIN_BOUNDS" ) {
			Size bin;
			Real lower;
			l >> tag >> bin >> tag >> lower;
			pair_bin_lower_bounds_[ bin ] = lower * lower;

		} else if ( tag == "PAIR" ) {
			Size bin;
			char na('?');
			l >> tag >> bin >> tag >> na >> tag >> aa >> tag >> enrich;
			Size const aa_index( aa_from_oneletter_code( aa ) );
			Size na_index( dna_aa_from_oneletter_code( na ) );
			if ( aa_index > num_canonical_aas || na_index < first_DNA_aa || na_index > last_DNA_aa ) {
				utility_exit_with_message( "bad line " + line );
			}
			na_index -= ( first_DNA_aa - 1 ); // now goes from 1 -> 4
			enrich = std::min( std::max( min_pair_enrichment_, enrich ), max_pair_enrichment_ );
			pair_stats_[ bin ][ na_index ][ aa_index ] = -1.0 * std::log( enrich );

		} else if ( tag == "MAX_NBRS" ) {
			l >> max_nbrs_;
			nbr_bin_.resize( max_nbrs_+1, 0 ); // vector0

		} else if ( tag == "NBR_BIN" ) {
			Size nbr_count(0);
			l >> tag >> nbr_count >> tag >> nbr_bin;
			nbr_bin_[ nbr_count ] = nbr_bin;
			if ( nbr_bin > env_stats_.size() ) {
				env_stats_.resize( nbr_bin, twenty_zeroes );
			}

		} else if ( tag == "ENV" ) {
			l >> tag >> nbr_bin >> tag >> aa >> tag >> enrich;
			Size const aa_index( aa_from_oneletter_code( aa ) );
			enrich = std::min( std::max( min_env_enrichment_, enrich ), max_env_enrichment_ );
			env_stats_[ nbr_bin ][ aa_index ] = -1.0 * std::log( enrich );

		} else {
			utility_exit_with_message( "unrecognized line "+line );

		}
		if ( l.fail() ) {
			utility_exit_with_message( "misparsed line "+line );
		}
	}
}



} // ns dna
} // ns scoring
} // ns core
