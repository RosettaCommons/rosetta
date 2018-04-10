// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/SimpleMetricDataLoaderCreator.hh
/// @brief  Creator classes for the default DataLoader classes, TaskOperationLoader and ScoreFunctionLoader
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_SimpleMetricLoaderCreator_hh
#define INCLUDED_protocols_jd2_parser_SimpleMetricLoaderCreator_hh

// Package headers
#include <protocols/parser/DataLoaderCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace parser {

class SimpleMetricLoaderCreator : public DataLoaderCreator
{
public:
	virtual DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
	virtual DerivedNameFunction schema_ct_naming_function() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace parser
} //namespace protocols

#endif
