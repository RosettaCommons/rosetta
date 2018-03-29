// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_filters/ShowConstraintsFilterCreator.hh
/// @brief  FilterCreators for the ShowConstraintsFilterCreator
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)


#ifndef INCLUDED_protocols_constraint_filters_ShowConstraintsFilterCreator_hh
#define INCLUDED_protocols_constraint_filters_ShowConstraintsFilterCreator_hh

// Package Headers
#include <protocols/filters/FilterCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace constraint_filters {

class ShowConstraintsFilterCreator : public protocols::filters::FilterCreator
{
public:
	virtual protocols::filters::FilterOP create_filter() const override;
	virtual std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd  ) const override;
};

}
}

#endif
