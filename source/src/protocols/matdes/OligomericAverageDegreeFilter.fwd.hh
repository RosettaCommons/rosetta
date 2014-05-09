// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/matdes/OligomericAverageDegreeFilter.fwd.hh
/// @brief
/// @author Neil King (neilking@u.washington.edu)


#ifndef INCLUDED_protocols_matdes_OligomericAverageDegreeFilter_fwd_hh
#define INCLUDED_protocols_matdes_OligomericAverageDegreeFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace matdes {

// Forward
class OligomericAverageDegreeFilter;

// Types
typedef utility::pointer::owning_ptr< OligomericAverageDegreeFilter >  OligomericAverageDegreeFilterOP;
typedef utility::pointer::owning_ptr< OligomericAverageDegreeFilter const >  OligomericAverageDegreeFilterCOP;

} // namespace matdes
} // namespace protocols

#endif
