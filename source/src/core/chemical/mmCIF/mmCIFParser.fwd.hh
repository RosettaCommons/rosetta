// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/chemical/mmCIF/mmCIFParser.fwd.hh
/// @author Steven Combs (steven.combs1@gmail.com)

#ifndef INCLUDED_core_chemical_mmCIF_mmCIFParser_fwd_hh
#define INCLUDED_core_chemical_mmCIF_mmCIFParser_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>
#include <cifparse/ISTable.h>

namespace core {
namespace chemical {
namespace mmCIF {

class mmCIFParser;

typedef utility::pointer::shared_ptr< CifFile > CifFileOP;
typedef utility::pointer::shared_ptr< CifParser > CifParserOP;

}
}
}

#endif

