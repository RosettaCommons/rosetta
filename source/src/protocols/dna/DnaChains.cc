// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaChains.cc
/// @author ashworth

#include <protocols/dna/DnaChains.hh>

#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <ObjexxFCL/format.hh> // I()

#include <utility/vector1.hh>


using utility::vector1;

namespace protocols {
namespace dna {

using namespace ObjexxFCL::format;
using namespace core;

DnaChains::DnaChains() : utility::pointer::ReferenceCount() {}

DnaChains::~DnaChains()= default;

DnaChains::DnaChains( DnaChains const & other )
: utility::pointer::ReferenceCount()
{
	positions_ = other.positions();
}

DnaChainsOP DnaChains::clone() const
{
	return DnaChainsOP( new DnaChains( *this ) );
}

bool
DnaChains::contains( Size index ) const
{
	if ( is_top(index) ) {
		runtime_assert( (*this)[ index ].top() == index );
		return true;
	}
	// bottom strand check
	 for ( auto const & position : positions_ ) {
		if ( index == position.second.bottom() ) return true;
	}
	return false;
}

void
DnaChains::print(
	pose::Pose const & pose,
	std::ostream & os
) const
{
	os << "There are " << positions_.size() << " dna positions:" << '\n';
	 for ( auto const & position : positions_ ) {

		Size const top_i( position.first );
		DnaPosition const & pos( position.second );
		runtime_assert( top_i == pos.top() );
		os << I( 4, top_i );
		if ( pose.pdb_info() ) {
			os << " (pdb " << pose.pdb_info()->chain( top_i ) << " "
				<< I( 4, pose.pdb_info()->number( top_i ) ) << ")";
		}
		os << pose.residue_type( top_i ).name3();
		if ( pos.paired() ) {
			Size const bot_i( pos.bottom() );
			os << "  <=>" << pose.residue_type( bot_i ).name3() << "  " << I( 4, bot_i );
			if ( pose.pdb_info() ) {
				os << " (pdb " << pose.pdb_info()->chain( bot_i ) << " "
					<< I( 4, pose.pdb_info()->number( bot_i ) ) << ")";
			}
		} else os << " (unpaired)";
		os << '\n';
	}
	os << std::endl;
}

} // namespace dna
} // namespace protocols
