// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/FileListResourceLocater.cc
/// @brief A resource locator that takes a space seperated list of files and concatenates them into a stringstream.  Useful mostly for concatenating pdbs
/// @author


//unit headers
#include <basic/resource_manager/locator/FileListResourceLocator.hh>
#include <basic/resource_manager/locator/FileListResourceLocatorCreator.hh>

#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>
#include <basic/resource_manager/locator/StringResourceStream.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

namespace basic {
namespace resource_manager {
namespace locator {

static Tracer file_list_tracer("basic.resource_manager.locator.FileListResourceLocator");

FileListResourceLocatorCreator::FileListResourceLocatorCreator() {}

FileListResourceLocatorCreator::~FileListResourceLocatorCreator(){}

ResourceLocatorOP
FileListResourceLocatorCreator::create_resource_locator() const
{
	return new FileListResourceLocator;
}

std::string
FileListResourceLocatorCreator::locator_type() const
{
	return "FileListResourceLocator";
}

FileListResourceLocator::FileListResourceLocator() : basic::resource_manager::ResourceLocator(),
	open_mode_()
{

}

FileListResourceLocator::~FileListResourceLocator()
{

}

FileListResourceLocator::FileListResourceLocator(
	FileListResourceLocator const & /*src*/
) : basic::resource_manager::ResourceLocator(),
	open_mode_()
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
	return "FileListResourceLocator";
}

ResourceStreamOP
FileListResourceLocator::locate_resource_stream(std::string const & locator_tag) const
{
	utility::vector1<std::string> path_vector(utility::string_split(locator_tag));

	StringResourceStreamOP string_stream = new StringResourceStream;

	for(utility::vector1<std::string>::const_iterator path_it = path_vector.begin();path_it != path_vector.end();++path_it)
	{
		FileStream new_file(*path_it);
		while(new_file.stream().good())
		{
			std::string line;
			getline(new_file.stream(),line);
			string_stream->fill(line + "\n");
		}
	}
	return string_stream;
}

void
FileListResourceLocator::parse_my_tag(utility::tag::TagCOP)
{

}

}
}
}
