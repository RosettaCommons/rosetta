// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictOperationsBase.cc
/// @brief  Base class for PoseMetricCalculator-using TaskOperations
/// @author Steven Lewis smlewi@gmail.com

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

// Utility Headers
//#include <core/types.hh>
#include <utility/vector1_bool.hh>
//#include <basic/Tracer.hh>

// C++ headers
#include <set>

#include <utility/vector1.hh>

//static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictOperationsBase" );

namespace protocols {
namespace toolbox {
namespace task_operations {

RestrictOperationsBase::RestrictOperationsBase() : parent() {}

RestrictOperationsBase::~RestrictOperationsBase() {}

void
RestrictOperationsBase::run_calculator(
	core::pose::Pose const & pose,
	std::string const & calculator,
	std::string const & calculation,
	utility::vector1_bool & residues ) const
{
	runtime_assert(residues.size() == pose.size());

	//find the set of residues
	typedef std::set< core::Size > SizeSet;
	basic::MetricValue< SizeSet > mv_sizeset;
	pose.metric(calculator, calculation, mv_sizeset);
	SizeSet const & sizeset(mv_sizeset.value());

	//insert this into the vector
	for ( SizeSet::const_iterator it(sizeset.begin()), end(sizeset.end()) ; it != end; ++it ) {
		//TR << *it << std::endl;
		residues[*it] = true;
	}

	return;
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
