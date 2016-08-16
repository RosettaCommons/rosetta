// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/fldsgn/filters/InterlockingAromaFilter.fwd.hh
/// @brief
/// @author Nobuyasu Koga (nobuyasu@uw.edu)


#ifndef INCLUDED_protocols_fldsgn_filters_InterlockingAromaFilter_fwd_hh
#define INCLUDED_protocols_fldsgn_filters_InterlockingAromaFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace fldsgn {
namespace filters {

// Forward
class InterlockingAromaFilter;

// Types
typedef utility::pointer::shared_ptr< InterlockingAromaFilter >  InterlockingAromaFilterOP;
typedef utility::pointer::shared_ptr< InterlockingAromaFilter const >  InterlockingAromaFilterCOP;

} // namespace protocols
} // namespace fldsgn
} // namespace filters

#endif
