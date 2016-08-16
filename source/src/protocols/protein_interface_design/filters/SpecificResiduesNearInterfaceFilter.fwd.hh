// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilter.fwd.hh
/// @brief  Forward header for SpecificResiduesNearInterface
/// @author Arpit Tandon
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilter_fwd_hh
#define INCLUDED_protocols_protein_interface_design_filters_SpecificResiduesNearInterfaceFilter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

// Forward
class SpecificResiduesNearInterfaceFilter;

// Types
typedef utility::pointer::shared_ptr< SpecificResiduesNearInterfaceFilter >  SpecificResiduesNearInterfaceFilterOP;
typedef utility::pointer::shared_ptr< SpecificResiduesNearInterfaceFilter const >  SpecificResiduesNearInterfaceFilterCOP;

} // namespace
} // namespace
} // namespace

#endif // Include guard
