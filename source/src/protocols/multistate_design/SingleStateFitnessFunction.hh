// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SingleStateFitnessFunction.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_SingleStateFitnessFunction_hh
#define INCLUDED_protocols_multistate_design_SingleStateFitnessFunction_hh

#include <protocols/multistate_design/SingleStateFitnessFunction.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

class SingleStateFitnessFunction : public utility::pointer::ReferenceCount {

public:
	SingleStateFitnessFunction() {}
	virtual ~SingleStateFitnessFunction() {}

	virtual core::Real evaluate( core::pose::Pose const & pose ) const;

private:
};

} // namespace multistate_design
} // namespace protocols

#endif
