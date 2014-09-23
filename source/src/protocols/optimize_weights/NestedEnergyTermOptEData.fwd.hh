// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file protocols/optimize_weights/NestedEnergyTermOptEData.fwd.hh
/// @brief Forward header for some special optE position data classes
/// @author Ron Jacak

#ifndef INCLUDED_protocols_optimize_weights_NestedEnergyTermOptEData_fwd_hh
#define INCLUDED_protocols_optimize_weights_NestedEnergyTermOptEData_fwd_hh


// utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace optimize_weights {

class NestedEnergyTermPNatAAOptEPositionData;
class NestedEnergyTermDDGMutationOptEData;

/// Position containers
typedef utility::pointer::shared_ptr< NestedEnergyTermPNatAAOptEPositionData > NestedEnergyTermPNatAAOptEPositionDataOP;
typedef utility::pointer::shared_ptr< NestedEnergyTermDDGMutationOptEData > NestedEnergyTermDDGMutationOptEDataOP;


} // namespace optimize_weights
} // namespace protocols


#endif
