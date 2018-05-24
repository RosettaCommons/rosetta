// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   RequirementCreator.hh
/// @brief  Forward Header for base class for RequirementCreators for the Requirement load-time factory registration scheme
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_requirements_SizeInSegmentsRequirementCreator_fwd_hh
#define INCLUDED_protocols_sewing_requirements_SizeInSegmentsRequirementCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace sewing  {
namespace requirements {

/// @brief Creates SizeInSegmentsRequirements for the AssemblyRequirementFactory
class SizeInSegmentsRequirementCreator;
typedef utility::pointer::shared_ptr< SizeInSegmentsRequirementCreator > SizeInSegmentsRequirementCreatorOP;
typedef utility::pointer::shared_ptr< SizeInSegmentsRequirementCreator const > SizeInSegmentsRequirementCreatorCOP;

} //namespace
} //namespace
} //namespace

#endif
