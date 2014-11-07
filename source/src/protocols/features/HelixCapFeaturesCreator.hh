// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/HelixCapFeaturesCreator.hh
/// @brief  
/// @author Ben Borgo
/// @author Jim Havranek

#ifndef INCLUDED_protocols_features_HelixCapFeaturesCreator_hh
#define INCLUDED_protocols_features_HelixCapFeaturesCreator_hh

// Unit Headers
#include <protocols/features/FeaturesReporterCreator.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace features {
	
/// @brief creator for the HelixCapFeaturesCreator class
class HelixCapFeaturesCreator : public FeaturesReporterCreator
{
public:
	HelixCapFeaturesCreator();
	virtual ~HelixCapFeaturesCreator();
	
	virtual FeaturesReporterOP create_features_reporter() const;
	virtual std::string type_name() const;
};
	
} //namespace features
} //namespace protocols

#endif
