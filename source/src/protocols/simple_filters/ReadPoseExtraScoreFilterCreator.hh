// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ReadPoseExtraScoreFilterCreator.hh
/// @brief  FilterCreators for the ReadPoseExtraScoreFilter
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_protocols_simple_filters_ReadPoseExtraScoreFilterCreator_hh
#define INCLUDED_protocols_simple_filters_ReadPoseExtraScoreFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>
#include <protocols/simple_filters/ReadPoseExtraScoreFilter.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace simple_filters {

class ReadPoseExtraScoreFilterCreator : public filters::FilterCreator
{
public:
	filters::FilterOP create_filter() const override {
		return protocols::filters::FilterOP( new ReadPoseExtraScoreFilter );
	}

	std::string keyname() const override {
		return ReadPoseExtraScoreFilter::class_name();
	}

	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override {
		return ReadPoseExtraScoreFilter::provide_xml_schema( xsd );
	}
};

}
}
#endif
