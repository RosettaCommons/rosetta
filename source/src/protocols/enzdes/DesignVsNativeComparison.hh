// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file stuff to compare designs against the native pdb
/// @brief
/// @author Florian Richter, floric@u.washington.edu


#ifndef INCLUDED_protocols_enzdes_DesignVsNativeComparison_hh
#define INCLUDED_protocols_enzdes_DesignVsNativeComparison_hh


// Unit headers

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentEnergy.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

//Utility Headers

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

class DesignVsNativeComparison;
typedef utility::pointer::shared_ptr< DesignVsNativeComparison > DesignVsNativeComparisonOP;


/// @brief class that holds a bunch of native poses and compares them
/// @brief against a given input pose on request
class DesignVsNativeComparison : public utility::pointer::ReferenceCount {

	//friend class EnzConstraintParameters;

public:

	DesignVsNativeComparison();

	~DesignVsNativeComparison() override;

	void
	compare_to_native(
		core::pose::Pose const & pose,
		utility::vector1< std::pair< std::string, std::string > > const & calculators,
		core::scoring::ScoreFunctionCOP scorefxn,
		utility::vector1< core::io::silent::SilentEnergy > & silent_Es
	);

private:

	std::map< std::string, core::pose::PoseOP > native_poses_;

}; //class DesignVsNativeComparison

} //protocols
} //enzdes


#endif
