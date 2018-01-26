// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/FragmentSetLoader.cc
/// @brief
/// @author

// unit headers
#include <core/fragment/FragmentSetLoader.hh>
#include <core/fragment/FragmentSetLoaderCreator.hh>

// package headers
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>

// Basic headers
#include <basic/resource_manager/loader_schemas.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

//C++ headers
#include <istream>

namespace core {
namespace fragment {

FragmentSetLoader::FragmentSetLoader() {}
FragmentSetLoader::~FragmentSetLoader() = default;

basic::resource_manager::ResourceCOP
FragmentSetLoader::create_resource(
	basic::resource_manager::ResourceManager &,
	utility::tag::TagCOP resource_tag,
	std::string const & input_id,
	std::istream & input_stream
) const
{
	Size top( resource_tag->getOption< Size >( "top", 0 ));
	Size ncopies( resource_tag->getOption< Size >( "ncopies", 1 ));
	bool annotate( resource_tag->getOption< bool >( "annotate", true ));

	FragmentIO fio( top, ncopies, annotate );

	FragSetOP frags = fio.read_data_from_stream( input_id, input_stream );
	return frags;
}

std::string
FragmentSetLoader::classname()
{
	return "FragmentSet";
}

void
FragmentSetLoader::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "top", xsct_non_negative_integer, "TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "ncopies", xsct_non_negative_integer, "TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "annotate", xsct_rosetta_bool, "TO DO", "true" );

	basic::resource_manager::resource_loader_xsd_type_definition_w_attributes( xsd,
		classname(), "Load a FragmentSet and hold it in the ResourceManager; the thread-safe way to load fragments"
		" a single time. All RosettaScript-able movers ought to switch to using the resource manager to load"
		" and eventually unload the fragments they rely on",
		attlist );
}

/// @details Return an owning pointer to a newly constructed default instance of FragmentSetLoader.
basic::resource_manager::ResourceLoaderOP FragmentSetLoaderCreator::create_resource_loader() const
{
	return basic::resource_manager::ResourceLoaderOP( new FragmentSetLoader() );
}

/// @details Return a string specifying the type of %ResourceLoader to create
std::string FragmentSetLoaderCreator::loader_type() const
{
	return FragmentSetLoader::classname();
}

void
FragmentSetLoaderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	FragmentSetLoader::provide_xml_schema( xsd );
}

} // namespace fragment
} // namespace core
