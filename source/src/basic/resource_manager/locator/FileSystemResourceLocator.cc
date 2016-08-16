// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/FileSystemResourceLocater.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocatorCreator.hh>

//project headers
#include <utility/tag/Tag.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>


//utility headers
#include <utility/vector1.hh>

//C++ headers
#include <istream>
#include <string>

namespace basic {
namespace resource_manager {
namespace locator {

using utility::tag::TagCOP;
using utility::file::FileName;
using utility::io::izstream;
using utility::vector1;
using std::string;
using std::endl;
using std::istream;
using basic::Tracer;

static Tracer TR("basic.resource_manager.locator.FileSystemResourceLocator");


///// FileSystemResourceLocatorCreator /////
FileSystemResourceLocatorCreator::FileSystemResourceLocatorCreator() {}

FileSystemResourceLocatorCreator::~FileSystemResourceLocatorCreator() {}

ResourceLocatorOP
FileSystemResourceLocatorCreator::create_resource_locator() const {
	return ResourceLocatorOP( new FileSystemResourceLocator );
}

string
FileSystemResourceLocatorCreator::locator_type() const {
	return "FileSystemResourceLocator";
}

///// FileStream //////

/// @detail This is private. The FileStream shouldn't be copied
FileStream::FileStream(
	FileStream const & src
) :
	ResourceStream( src ),
	stream_()
{}

/// @detail If you use this constructor be sure to use the open
///function before accessing the stream
FileStream::FileStream() :
	ResourceStream(),
	stream_()
{}

FileStream::FileStream(
	string const & filename,
	std::ios_base::openmode open_mode
) :
	stream_(filename, open_mode)
{
	if ( !stream_ ) {
		vector1<string> alternative_search_paths(
			izstream::get_alternative_search_paths());

		TR << "Unable to open file '" << filename << "' at any of the following paths:" << endl;
		TR << "\t" << filename << endl;
		for (
				vector1<string>::const_iterator
				p=alternative_search_paths.begin(), pe=alternative_search_paths.end();
				p != pe; ++p ) {
			TR << "\t" << *p << platform::file::PATH_SEPARATOR << filename << endl;
			if ( utility::file::file_extension( filename) != "gz" ) {
				TR << "\t" << *p << filename << ".gz" << endl;
			}
		}
		throw utility::excn::EXCN_FileNotFound(filename);
	}
}

FileStream::~FileStream() {}

void
FileStream::open(
	string const & filename,
	std::ios_base::openmode open_mode
) {
	stream_.open(filename, open_mode);
	if ( !stream_ ) {
		throw utility::excn::EXCN_FileNotFound(filename);
	}

}

istream &
FileStream::stream() {
	return stream_;
}


///// FileSystemResourceLocator /////

FileSystemResourceLocator::FileSystemResourceLocator(
	std::ios_base::openmode open_mode
) :
	basic::resource_manager::ResourceLocator(),
	open_mode_(open_mode),
	base_path_("")
{}


FileSystemResourceLocator::FileSystemResourceLocator(
	FileSystemResourceLocator const & src
) :
	basic::resource_manager::ResourceLocator(),
	open_mode_(src.open_mode_),
	base_path_(src.base_path_)
{}

void
FileSystemResourceLocator::show(
	std::ostream & out
) const {
	out
		<< "FileSystemResourceLocator:" << endl
		<< "  open_mode:"
		<< (std::ios_base::app && open_mode_ ? " append" : "")
		<< (std::ios_base::ate && open_mode_ ? " at_end" : "")
		<< (std::ios_base::binary && open_mode_ ? " binary" : "")
		<< (std::ios_base::in && open_mode_ ? " input" : "")
		<< (std::ios_base::out && open_mode_ ? " output" : "")
		<< (std::ios_base::trunc && open_mode_ ? " truncate" : "")
		<< "  base search path for resources: "
		<< (base_path_ == "" ? "none" : base_path_)
		<< std::endl;
}

//std::ostream &
//FileSystemResourceLocator::operator<< (
// std::ostream & out,
// const FileSystemResourceLocator & file_system_resource_locator
//) {
// file_system_resource_locator.show(out);
// return out;
//}

std::string
FileSystemResourceLocator::type() const
{
	return "FileSystemResourceLocator";
}

void
FileSystemResourceLocator::set_open_mode(
	std::ios_base::openmode open_mode
)
{
	open_mode_ = open_mode;
}

std::ios_base::openmode
FileSystemResourceLocator::get_open_mode() const {
	return open_mode_;
}

FileSystemResourceLocator::~FileSystemResourceLocator() {}

/// @brief
ResourceStreamOP
FileSystemResourceLocator::locate_resource_stream(
	string const & locator_tag
) const {
	// Concatenate base_path_ and the locator tag to generate the appropriate filename.
	std::stringstream fully_specified_locator_tag;
	fully_specified_locator_tag << base_path_ << locator_tag;
	return ResourceStreamOP( new FileStream( fully_specified_locator_tag.str(), open_mode_ ) );
}

/// @details Set the value for base_path if specified in the ResourceDefintionFile.
void
FileSystemResourceLocator::parse_my_tag(
	TagCOP tag
)
{
	if ( tag && tag->hasOption("base_path") ) {
		std::stringstream base_path_with_trailing_space;
		base_path_with_trailing_space << tag->getOption<string>("base_path") << "/";
		base_path_ = base_path_with_trailing_space.str();
	}
}

} // namespace locator
} // namespace resource_manager
} // namespace basic
