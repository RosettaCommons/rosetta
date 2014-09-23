// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DatabaseEntryWorkUnit.fwd.hh
///
/// @brief A work unit that runs a database query, processes the results, and returns a string (presumably a database insert statement)

/// @author Tim Jacobs

#ifndef INCLUDED_protocols_wum_DatabaseEntryWorkUnit_fwd_hh
#define INCLUDED_protocols_wum_DatabaseEntryWorkUnit_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace wum{

class DatabaseEntryWorkUnit;
typedef utility::pointer::shared_ptr< DatabaseEntryWorkUnit > DatabaseEntryWorkUnitOP;
typedef utility::pointer::shared_ptr< DatabaseEntryWorkUnit const > DatabaseEntryWorkUnitCOP;

}//namespace wum
}//namespace protocols

#endif
