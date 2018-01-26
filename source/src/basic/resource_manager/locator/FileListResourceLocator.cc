// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/FileListResourceLocater.cc
/// @brief A resource locator that takes a space seperated list of files and concatenates them into a stringstream.  Useful mostly for concatenating pdbs
/// @author


//unit headers
#include <basic/resource_manager/locator/FileListResourceLocator.hh>
#include <basic/resource_manager/locator/FileListResourceLocatorCreator.hh>

// Package headers
#include <basic/resource_manager/locator/locator_schemas.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>
#include <basic/resource_manager/locator/StringResourceStream.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace basic {
namespace resource_manager {
namespace locator {

static Tracer file_list_tracer("basic.resource_manager.locator.FileListResourceLocator");

FileListResourceLocatorCreator::FileListResourceLocatorCreator() = default;

FileListResourceLocatorCreator::~FileListResourceLocatorCreator()= default;

ResourceLocatorOP
FileListResourceLocatorCreator::create_resource_locator() const
{
	return ResourceLocatorOP( new FileListResourceLocator );
}

std::string
FileListResourceLocatorCreator::locator_type() const
{
	return FileListResourceLocator::classname();
}

void
FileListResourceLocatorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	return FileListResourceLocator::provide_xml_schema( xsd );
}

FileListResourceLocator::FileListResourceLocator() :
	basic::resource_manager::ResourceLocator(),
	open_mode_( std::ios_base::in )
{

}

FileListResourceLocator::~FileListResourceLocator() = default;

FileListResourceLocator::FileListResourceLocator(
	FileListResourceLocator const & /*src*/
) :
	basic::resource_manager::ResourceLocator(),
	open_mode_( std::ios_base::in )
{

}

void
FileListResourceLocator::show(std::ostream & out) const
{
	out << "FileListResourceLocator: " <<std::endl;
}

std::string
FileListResourceLocator::type() const
{
	return classname();
}

std::string
FileListResourceLocator::classname()
{
	return "FileListResourceLocator";
}

void
FileListResourceLocator::set_open_mode(
	std::ios_base::openmode open_mode
)
{
	open_mode_ = open_mode;
}

std::ios_base::openmode
FileListResourceLocator::get_open_mode() const {
	return open_mode_;
}

ResourceStreamOP
FileListResourceLocator::locate_resource_stream(std::string const & input_id) const
{
	utility::vector1<std::string> path_vector(utility::string_split(input_id));

	StringResourceStreamOP string_stream( new StringResourceStream );

	for ( utility::vector1<std::string>::const_iterator
			path_it = path_vector.begin();
			path_it != path_vector.end();
			++path_it ) {
		FileStream new_file( *path_it, open_mode_ );
		while ( new_file.stream().good() ) {
			std::string line;
			getline(new_file.stream(),line);
			string_stream->fill(line + "\n");
		}
	}
	return string_stream;
}

void
FileListResourceLocator::parse_my_tag(utility::tag::TagCOP) {}

void
FileListResourceLocator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
)
{
	using namespace utility::tag;
	AttributeList attrs;
	xsd_type_definition_w_attributes( xsd, classname(), "The file list resource locator will"
		" construct a single stream from the contents of one or more files listed in the 'input_id'"
		" of the Resource that needs to be constructed. These files should be separated by spaces."
		" This could be useful for constructing a single Pose from two separate chains, e.g.. It"
		" will not search for these files beyond the paths given in the input_id.", attrs );
}

}
}
}
