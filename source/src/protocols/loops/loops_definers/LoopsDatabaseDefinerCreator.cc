// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDatabaseDefinerCreator.cc
/// @brief  LoopsDatabaseDefinerCreator for the LoosDatabaseDefiner load-time factory registration scheme
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDatabaseDefinerCreator.hh>
#include <protocols/loops/loops_definers/LoopsDatabaseDefiner.hh>

using std::string;

namespace protocols {
namespace loops {
namespace loops_definers {

LoopsDatabaseDefinerCreator::LoopsDatabaseDefinerCreator() {}
LoopsDatabaseDefinerCreator::~LoopsDatabaseDefinerCreator() {}

LoopsDefinerOP
LoopsDatabaseDefinerCreator::create_loops_definer() const {
	return LoopsDefinerOP( new LoopsDatabaseDefiner );
}

string
LoopsDatabaseDefinerCreator::type_name() const {
	// return "LoopsDatabase";
	return LoopsDatabaseDefiner::class_name();
}

void LoopsDatabaseDefinerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopsDatabaseDefiner::provide_xml_schema( xsd );
}

} //namespace
} //namespace
} //namespace
