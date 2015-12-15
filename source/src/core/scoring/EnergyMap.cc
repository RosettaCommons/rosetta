// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/EnergyMap.cc
/// @brief  Vector of scores implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit Headers
#include <core/scoring/EnergyMap.hh>

// Package Headers
#include <core/scoring/ScoreType.hh>

#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION


namespace core {
namespace scoring {

void
EMapVector::print() const
{
	std::cout << "( ";
	for ( int ii = 0; ii < n_score_types; ++ii ) {
		if ( ii != 0 ) std::cout << ", ";
		std::cout << name_from_score_type( ScoreType( ii + 1 ) ) << "=" << map_[ ii ];
		if ( ii % 3 == 2 && ii != n_score_types ) std::cout << std::endl;
	}
	std::cout << ")" << std::endl;
}


void
EMapVector::show_if_nonzero_weight( std::ostream & out, EMapVector const & weights ) const
{
	for ( int ii = 1; ii <= n_score_types; ++ii ) {
		if ( weights[ ScoreType(ii) ] != 0.0 ) {
			Real const val( operator[]( ScoreType(ii) ) );
			out << ' ' << ScoreType(ii) << ": " << ObjexxFCL::format::F(9,3,val);
		}
	}
}


void
EMapVector::show_weighted( std::ostream & out, EMapVector const & weights ) const
{
	for ( int ii = 1; ii <= n_score_types; ++ii ) {
		Real const weight( weights[ ScoreType(ii) ] );
		if ( weight  != 0.0 ) {
			Real const val( operator[]( ScoreType(ii) ) );
			out << ' ' << ScoreType(ii) << ": " << ObjexxFCL::format::F(9,3, weight * val );
		}
	}
}


std::string
EMapVector::weighted_string_of( EMapVector const & weights ) const
{
	std::ostringstream os;
	show_weighted( os, weights );
	return os.str();
}


} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::EMapVector::save( Archive & arc ) const {
	arc( CEREAL_NVP( map_ ) ); // Real [361]
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::EMapVector::load( Archive & arc ) {
	arc( map_ ); // Real [361]
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::EMapVector );

#endif // SERIALIZATION
