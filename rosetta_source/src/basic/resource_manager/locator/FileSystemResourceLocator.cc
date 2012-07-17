// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/FileSystemResourceLocater.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/locator/FileSystemResourceLocator.hh>
#include <basic/resource_manager/locator/FileSystemResourceLocatorCreator.hh>

//project headers
#include <utility/tag/Tag.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>


//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <istream>
#include <string>

namespace basic {
namespace resource_manager {
namespace locator {

using utility::tag::TagPtr;
using utility::file::FileName;
using utility::io::izstream;
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
	return new FileSystemResourceLocator;
}

string
FileSystemResourceLocatorCreator::locator_type() const {
	return "FileSystemResourceLocator";
}

///// FileStream //////

///@detail This is private. The FileStream shouldn't be copied
FileStream::FileStream(
	FileStream const & src
) :
	ResourceStream( src ),
	stream_()
{}

///@detail If you use this constructor be sure to use the open
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
{}

FileStream::~FileStream() {}

void
FileStream::open(
	string const & filename,
	std::ios_base::openmode open_mode
) {
	stream_.open(filename, open_mode);
}

istream &
FileStream::stream() {
	return stream_;
}


///// FileSystemResourceLocator /////

FileSystemResourceLocator::FileSystemResourceLocator(
	std::ios_base::openmode open_mode
) :
	open_mode_(open_mode)
{}


FileSystemResourceLocator::FileSystemResourceLocator(
  FileSystemResourceLocator const & src
) :
	open_mode_(src.open_mode_)
{}

void
FileSystemResourceLocator::set_open_mode(
	std::ios_base::openmode open_mode
) {
	open_mode_ = open_mode;
}

std::ios_base::openmode
FileSystemResourceLocator::get_open_mode() const {
	return open_mode_;
}

FileSystemResourceLocator::~FileSystemResourceLocator() {}

/// @brief Create a ResourceStream object from the given resource
/// source, so that its stream can be passed to the ResourceLoader
ResourceStreamOP
FileSystemResourceLocator::locate_resource_stream(
	string const & locator_tag
) const {
	return new FileStream( locator_tag );
}

void
FileSystemResourceLocator::parse_my_tag(
	TagPtr
)
{}

} // namespace locator
} // namespace resource_manager
} // namespace basic
