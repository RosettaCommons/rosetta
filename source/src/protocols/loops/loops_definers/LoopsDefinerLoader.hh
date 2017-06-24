// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsDefinerLoader.hh
/// @brief  Declartion of the XML loops_definers's LoopsDefinerLoader class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsDefinerLoader_hh
#define INCLUDED_protocols_loops_loops_definers_LoopsDefinerLoader_hh

// Package Headers
#include <protocols/parser/DataLoader.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace loops {
namespace loops_definers {

/// @brief A class for loading Loops data into the XML parser's basic::datacache::DataMap.
class LoopsDefinerLoader : public parser::DataLoader
{
public:
	LoopsDefinerLoader();
	virtual ~LoopsDefinerLoader();

	/// @brief The LoopsDefinerLoader will load named task operations into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

	static std::string loader_name();
	static std::string loop_def_loader_ct_namer( std::string const & element_name );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
};

} //namespace
} //namespace
} //namespace

#endif
