// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/optimize_weights/DDGBindOptEData.fwd.hh
/// @brief  Forward header for a optE position data class that holds interface ddG information
/// @author Ron Jacak


#ifndef INCLUDED_protocols_optimize_weights_DDGBindOptEData_fwd_hh
#define INCLUDED_protocols_optimize_weights_DDGBindOptEData_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace optimize_weights {

class DDGBindOptEData;
typedef utility::pointer::shared_ptr< DDGBindOptEData > DDGBindOptEDataOP;


} // namespace optimize_weights
} // namespace protocols

#endif // INCLUDED_protocols_optimize_weights_DDGBindOptEData_FWD_HH
