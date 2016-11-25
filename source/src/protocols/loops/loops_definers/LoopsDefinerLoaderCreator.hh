// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDefinerLoaderCreator.hh
/// @brief  Creator classe for the LoopsDefinerLoader classe
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsDefinerLoaderCreator_hh
#define INCLUDED_protocols_loops_loops_definers_LoopsDefinerLoaderCreator_hh

// Package headers
#include <protocols/jd2/parser/DataLoaderCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace loops {
namespace loops_definers {


class LoopsDefinerLoaderCreator : public jd2::parser::DataLoaderCreator
{
public:
	virtual jd2::parser::DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
	virtual DerivedNameFunction schema_ct_naming_function() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};


} //namespace
} //namespace
} //namespace

#endif
