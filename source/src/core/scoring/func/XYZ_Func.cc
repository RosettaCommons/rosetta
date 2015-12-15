// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

/// Unit headers
#include <core/scoring/func/XYZ_Func.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace func {

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/// helper classes to reuse code
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


XYZ_Func::~XYZ_Func(){}

ResidueXYZ::ResidueXYZ( Residue const & rsd_in ): rsd_( rsd_in ) {}

Vector const &
ResidueXYZ::operator()( AtomID const & id ) const
{
	debug_assert( id.rsd() == rsd_.seqpos() && id.atomno() <= rsd_.natoms() );
	return rsd_.xyz( id.atomno() );
}


conformation::Residue const &
ResidueXYZ::residue( Size ASSERT_ONLY( seqpos ) ) const {
	debug_assert( (Size) rsd_.seqpos() == seqpos );
	return rsd_;
}


ResidueXYZ::~ResidueXYZ(){}

ResiduePairXYZ::ResiduePairXYZ( Residue const & rsd1_in, Residue const & rsd2_in )
:
	rsd1_( rsd1_in ),
	rsd2_( rsd2_in )
{}

Vector const &
ResiduePairXYZ::operator()( AtomID const & id ) const
{
	Residue const & rsd( ( id.rsd() == rsd1_.seqpos() ) ? rsd1_ : rsd2_ );
	debug_assert( id.rsd() == rsd.seqpos() && id.atomno() <= rsd.natoms() );
	return rsd.xyz( id.atomno() );
}


conformation::Residue const &
ResiduePairXYZ::residue( Size seqpos ) const {
	debug_assert( (Size) rsd1_.seqpos() == seqpos || (Size) rsd2_.seqpos() == seqpos );
	return (Size) rsd1_.seqpos() == seqpos ? rsd1_ : rsd2_;
}

ResiduePairXYZ::~ResiduePairXYZ() {}

ConformationXYZ::ConformationXYZ( Conformation const & conformation_in ): conformation_( conformation_in ) {}


Vector const &
ConformationXYZ::operator()( AtomID const & id ) const
{
	return conformation_.xyz( id );
}


conformation::Residue const &
ConformationXYZ::residue( Size seqpos ) const {
	return conformation_.residue( seqpos );
}

ConformationXYZ::~ConformationXYZ() {}

} // constraints
} // scoring
} // core
