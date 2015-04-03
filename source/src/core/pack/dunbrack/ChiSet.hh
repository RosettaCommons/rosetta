// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/ChiSet.hh
/// @brief  Typedefs and forward declarations for class DunbrackRotamer
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_pack_dunbrack_ChiSet_hh
#define INCLUDED_core_pack_dunbrack_ChiSet_hh

// Package headers
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh> // where ChiVector and RotVector live

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace pack {
namespace dunbrack {

class ChiSet : public utility::pointer::ReferenceCount
{
public:
	ChiSet():
		chi( /* 0 */ ),
		rot( /* 0 */ ),
		ex_chi_steps( /* 0 */ ),
		probability( 0.0 )
	{}

	ChiSet( Size nchi ) :
		chi( nchi, 0 ),
		rot( nchi, 0 ),
		ex_chi_steps( nchi, 0 ),
		probability( 0.0 )
	{}

	ChiSet( ChiSet const & chi_set ) :
		ReferenceCount(),
		chi( chi_set.chi ),
		rot( chi_set.rot ),
		ex_chi_steps( chi_set.ex_chi_steps ),
		probability( chi_set.probability )
	{}

	ChiVector chi;
	RotVector rot;
	utility::vector1< Real > ex_chi_steps;

	/// total rperc for this chiset, including
	Real probability;

};


}
}
}

#endif //INCLUDED_core_pack_dunbrack_ChiSet_FWD_HH


