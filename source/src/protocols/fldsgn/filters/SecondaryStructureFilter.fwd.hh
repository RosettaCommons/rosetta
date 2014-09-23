// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/fldsgn/filters/SecondaryStructureFilter.fwd.hh
/// @brief
/// @author Nobuyasu Koga (nobuyasu@uw.edu)


#ifndef INCLUDED_protocols_fldsgn_filters_SecondaryStructureFilter_fwd_hh
#define INCLUDED_protocols_fldsgn_filters_SecondaryStructureFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace fldsgn {
namespace filters {

// Forward
class SecondaryStructureFilter;

// Types
typedef utility::pointer::shared_ptr< SecondaryStructureFilter >  SecondaryStructureFilterOP;
typedef utility::pointer::shared_ptr< SecondaryStructureFilter const >  SecondaryStructureFilterCOP;

} // namespace protocols
} // namespace fldsgn
} // namespace filters

#endif
