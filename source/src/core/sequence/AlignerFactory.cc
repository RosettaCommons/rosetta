// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/sequence/AlignerFactory.cc
/// @author James Thompson

// Unit headers
#include <core/sequence/Aligner.hh>
#include <core/sequence/AlignerFactory.hh>

// Package headers
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/MCAligner.hh>

#include <utility/exit.hh>

using namespace core::sequence;

namespace core {
namespace sequence {

AlignerOP
AlignerFactory::get_aligner(
	std::string const & type
) {
	if ( type == "local" ) {
		return AlignerOP( new SWAligner );
	} else if ( type == "global" ) {
		return AlignerOP( new NWAligner );
	} else if ( type == "mc" ) {
		return AlignerOP( new MCAligner );
	} else {
		utility_exit_with_message(
			"AlignerFactory: unknown Aligner type: " + type
		);
		return NULL;
	}
}

} // namespace sequence
} // namespace core
