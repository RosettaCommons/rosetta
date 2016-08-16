// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/CalculatorFilter.fwd.hh
/// @brief  forward declaration for CalculatorFilter
/// @author Rocco Moretti (rmoretti@u.washington.edu)


#ifndef INCLUDED_protocols_filters_CalculatorFilter_fwd_hh
#define INCLUDED_protocols_filters_CalculatorFilter_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace filters {

// Forward
class CalculatorFilter;

// Types
typedef utility::pointer::shared_ptr< CalculatorFilter >  CalculatorFilterOP;
typedef utility::pointer::shared_ptr< CalculatorFilter const >  CalculatorFilterCOP;

} // namespace filters
} // namespace protocols

#endif
