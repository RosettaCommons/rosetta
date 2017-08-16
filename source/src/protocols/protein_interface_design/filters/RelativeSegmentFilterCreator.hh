// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_interface_design/filters/RelativeSegmentFilterCreator.hh
/// @brief  FilterCreators for the filters defined in RelativeSegmentFilterCreator.hh
/// @author Sarel

#ifndef INCLUDED_protocols_protein_interface_design_filters_RelativeSegmentFilterCreator_hh
#define INCLUDED_protocols_protein_interface_design_filters_RelativeSegmentFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace protein_interface_design {
namespace filters {

class RelativeSegmentFilterCreator : public protocols::filters::FilterCreator
{
public:
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
} //namespace protein_interface_design
} //namespace protocols

#endif
