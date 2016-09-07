// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pose_reporters/BasicPosePropertyReporterCreators.hh
/// @brief  Collection of PosePropertySelector creators for basic pose reporters
/// @author Luki Goldschmidt <luki@uw.edu>

#ifndef INCLUDED_protocols_pose_reporters_BasicPosePropertyReporterCreators_hh
#define INCLUDED_protocols_pose_reporters_BasicPosePropertyReporterCreators_hh

#include <protocols/rosetta_scripts/PosePropertyReporterCreator.hh>
#include <string>

namespace protocols {
namespace pose_reporters {

class EnergyReporterCreator : public protocols::rosetta_scripts::PosePropertyReporterCreator {
public:
	protocols::rosetta_scripts::PosePropertyReporterOP create_reporter() const override;
	std::string keyname() const override { return "EnergyReporter"; }
};

class FilterReporterCreator : public protocols::rosetta_scripts::PosePropertyReporterCreator {
public:
	protocols::rosetta_scripts::PosePropertyReporterOP create_reporter() const override;
	std::string keyname() const override { return "FilterReporter"; }
};

class RMSDReporterCreator : public protocols::rosetta_scripts::PosePropertyReporterCreator {
public:
	protocols::rosetta_scripts::PosePropertyReporterOP create_reporter() const override;
	std::string keyname() const override { return "RMSDReporter"; }
};

}
}

#endif

