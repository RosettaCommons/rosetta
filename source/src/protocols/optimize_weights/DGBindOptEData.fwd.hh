// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/optimize_weights/DGBindOptEData.fwd.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_optimize_weights_DGBindOptEData_fwd_hh
#define INCLUDED_protocols_optimize_weights_DGBindOptEData_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace optimize_weights {


class DGBindOptEData; // fwd declaration
typedef utility::pointer::shared_ptr< DGBindOptEData > DGBindOptEDataOP;
typedef utility::pointer::shared_ptr< DGBindOptEData const > DGBindOptEDataCOP;


} // namespace optimize_weights
} // namespace protocols

#endif // INCLUDED_protocols_optimize_weights_DGBindOptEData_FWD_HH
