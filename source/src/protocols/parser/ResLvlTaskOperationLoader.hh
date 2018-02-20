// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/parser/ResLvlTaskOperationLoader.hh
/// @brief  Declartion of the XML parser's ResLvlTaskOperationLoader class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_ResLvlTaskOperationLoader_hh
#define INCLUDED_protocols_jd2_parser_ResLvlTaskOperationLoader_hh

// Package Headers
#include <protocols/parser/DataLoader.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace parser {

/// @brief A class for loading arbitrary data into the XML parser's basic::datacache::DataMap.
class ResLvlTaskOperationLoader : public DataLoader
{
public:
	ResLvlTaskOperationLoader();
	virtual ~ResLvlTaskOperationLoader();

	/// @brief The ResLvlTaskOperationLoader will load named residue-level-task operations
	/// into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

	static std::string loader_name();
	static std::string res_lvl_task_op_loader_ct_namer( std::string const & element_name );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

};

} //namespace parser
} //namespace protocols

#endif
