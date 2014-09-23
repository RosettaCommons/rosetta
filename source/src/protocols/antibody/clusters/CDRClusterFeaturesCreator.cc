
#include "CDRClusterFeaturesCreator.hh"

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterFeaturesCreator.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/CDRClusterFeaturesCreator.hh>
#include <protocols/antibody/clusters/CDRClusterFeatures.hh>


namespace protocols {
namespace antibody {
namespace clusters {
	using namespace protocols::features;
	
CDRClusterFeaturesCreator::CDRClusterFeaturesCreator(){}
CDRClusterFeaturesCreator::~CDRClusterFeaturesCreator(){}

FeaturesReporterOP
CDRClusterFeaturesCreator::create_features_reporter() const {
	return FeaturesReporterOP( new CDRClusterFeatures );
}

std::string
CDRClusterFeaturesCreator::type_name() const {
	return "CDRClusterFeatures";
}
		
	
}
}
}
