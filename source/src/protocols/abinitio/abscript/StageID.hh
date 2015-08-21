// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/StageID.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_StageID_hh
#define INCLUDED_protocols_abinitio_abscript_StageID_hh

// Unit Headers

// Package headers

// Project headers

// Utility Headers
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

// C++ Headers

namespace protocols {
namespace abinitio {
namespace abscript {

enum StageID {
	I = 1,
	II,
	IIIa, IIIb,
	IVa, IVb,
	END
};

typedef utility::vector1< StageID > StageIDs;

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_StageID_hh
