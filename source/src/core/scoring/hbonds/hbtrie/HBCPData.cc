// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/CountPairData_1_1.cc
/// @brief  CountPair Data for a residue with one connection point
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/hbtrie/HBCPData.hh>

// STL Headers
#include <iostream>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

HBCPData::HBCPData() : avoid_sc_hbonds_( false ), is_sc_( false ) {}

void
HBCPData::print( std::ostream & os ) const
{
	os << "HBCPData: avoid_sc_hbonds_= " << avoid_sc_hbonds_ << "; is_sc_= " << is_sc_ << ";";
}

//void
//HBCPData::set_count_pair_data_to_use(
// Size /*connection_id*/
//) const
//{}

std::ostream & operator << ( std::ostream & os, HBCPData const & hbcpdat )
{
	hbcpdat.print( os );
	return os;
}

} // namespace hbtrie
} // namespace hbonds
} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::hbtrie::HBCPData::save( Archive & arc ) const {
	arc( CEREAL_NVP( avoid_sc_hbonds_ ) ); // _Bool
	arc( CEREAL_NVP( is_sc_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::hbtrie::HBCPData::load( Archive & arc ) {
	arc( avoid_sc_hbonds_ ); // _Bool
	arc( is_sc_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::hbtrie::HBCPData );
#endif // SERIALIZATION
