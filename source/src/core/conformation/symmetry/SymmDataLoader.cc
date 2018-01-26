// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmDataLoader.cc
/// @brief  load the SymmData data-structure, which is used to configure symmetric poses.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

//unit headers
#include <core/conformation/symmetry/SymmDataLoader.hh>
#include <core/conformation/symmetry/SymmDataLoaderCreator.hh>
#include <core/conformation/symmetry/SymmData.hh>

//package headers
#include <basic/resource_manager/types.hh>
#include <basic/resource_manager/loader_schemas.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//C++ headers
#include <istream>

namespace core {
namespace conformation {
namespace symmetry {

SymmDataLoader::SymmDataLoader() {}
SymmDataLoader::~SymmDataLoader() = default;

basic::resource_manager::ResourceCOP
SymmDataLoader::create_resource(
	basic::resource_manager::ResourceManager &,
	utility::tag::TagCOP,
	std::string const &,
	std::istream & symstream
) const {
	SymmDataOP symm_data( new SymmData() );
	symm_data->read_symmetry_data_from_stream( symstream );
	return symm_data;
}

std::string
SymmDataLoader::classname()
{
	return "SymmData";
}

void
SymmDataLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd,
		classname(), "Load a SymmData object in from the indicated input_id",
		attlist );
}

//// SymmDataLoaderCreator
basic::resource_manager::ResourceLoaderOP
SymmDataLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new SymmDataLoader() );
}

std::string SymmDataLoaderCreator::loader_type() const
{
	return SymmDataLoader::classname();
}

void
SymmDataLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SymmDataLoader::provide_xml_schema( xsd );
}

} // namespace
} // namespace
} // namespace
