// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/architects/util.hh
/// @brief utility functions for architects
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_util_hh
#define INCLUDED_protocols_denovo_design_architects_util_hh

// Protocol headers
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace architects {

DeNovoArchitectCOP
retrieve_denovo_architect( std::string const & architect_name, basic::datacache::DataMap & data );

void
store_denovo_architect( DeNovoArchitectOP architect, basic::datacache::DataMap & data );

} //protocols
} //denovo_design
} //architects


#endif	//protocols/denovo_design/architects_util_hh
