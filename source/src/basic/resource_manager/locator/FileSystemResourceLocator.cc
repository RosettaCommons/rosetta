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

// Package headers
#include <basic/resource_manager/locator/locator_schemas.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

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
FileSystemResourceLocatorCreator::FileSystemResourceLocatorCreator() = default;

FileSystemResourceLocatorCreator::~FileSystemResourceLocatorCreator() = default;

ResourceLocatorOP
FileSystemResourceLocatorCreator::create_resource_locator() const {
	return ResourceLocatorOP( new FileSystemResourceLocator );
}

string
FileSystemResourceLocatorCreator::locator_type() const {
	return FileSystemResourceLocator::classname();
}

void
FileSystemResourceLocatorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FileSystemResourceLocator::provide_xml_schema( xsd );
}

///// FileStream //////

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
		throw CREATE_EXCEPTION(utility::excn::FileNotFound, filename);
	}
}

FileStream::~FileStream() = default;

void
FileStream::open(
	string const & filename,
	std::ios_base::openmode open_mode
) {
	stream_.open(filename, open_mode);
	if ( !stream_ ) {
		throw CREATE_EXCEPTION(utility::excn::FileNotFound, filename);
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
	search_paths_( 1, "./" )
{
	if ( basic::options::option[ basic::options::OptionKeys::in::path::path ].user() ) {
		for ( std::string path : basic::options::option[ basic::options::OptionKeys::in::path::path ]() ) {
			if ( path.size() > 0 && path[ path.size()-1 ] == '/' ) {
				search_paths_.push_back( path );
			} else {
				search_paths_.push_back( path + "/" );
			}
		}
	}
}


FileSystemResourceLocator::FileSystemResourceLocator(
	FileSystemResourceLocator const & src
) :
	basic::resource_manager::ResourceLocator(),
	open_mode_(src.open_mode_),
	search_paths_( 1, "./" )
{}

FileSystemResourceLocator::~FileSystemResourceLocator() = default;

void
FileSystemResourceLocator::set_search_paths(
	utility::vector1< std::string > const & search_paths
)
{
	search_paths_ = search_paths;
}

void
FileSystemResourceLocator::show(
	std::ostream & out
) const {
	out
		<< "FileSystemResourceLocator:" << endl
		<< "  open_mode:"
		<< ((std::ios_base::app & open_mode_) ? " append" : "")
		<< ((std::ios_base::ate & open_mode_) ? " at_end" : "")
		<< ((std::ios_base::binary & open_mode_) ? " binary" : "")
		<< ((std::ios_base::in & open_mode_) ? " input" : "")
		<< ((std::ios_base::out & open_mode_) ? " output" : "")
		<< ((std::ios_base::trunc & open_mode_) ? " truncate" : "")
		<< "  search paths for resources:";
	for ( auto const & search_path : search_paths_ ) {
		out << " " << search_path;
	}
	out << std::endl;
}

std::string
FileSystemResourceLocator::type() const
{
	return classname();
}

std::string
FileSystemResourceLocator::classname()
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

/// @brief
ResourceStreamOP
FileSystemResourceLocator::locate_resource_stream(
	string const & input_id
) const {
	// Concatenate base_path_ and the locator tag to generate the appropriate filename.
	for ( std::string const & search_path : search_paths_ ) {
		std::string fname = search_path + input_id;
		std::ifstream ifstr( fname.c_str() );
		if ( ! ifstr.fail() ) {
			return ResourceStreamOP( new FileStream( fname, open_mode_ ) );
		}
	}

	std::ostringstream oss;
	oss << "Error in trying to locate file '" << input_id << "' in any of the following locations:\n";
	for ( std::string const & search_path : search_paths_ ) {
		oss << "    " << search_path + input_id << "\n";
	}
	oss << "Note that it is possible that the file exists by that you do not have permission to read it.\n";
	throw CREATE_EXCEPTION( utility::excn::Exception, oss.str() );

	// appease compiler
	return ResourceStreamOP();
}

/// @details Set the value for base_path if specified in the ResourceDefinitionFile.
void
FileSystemResourceLocator::parse_my_tag(
	TagCOP tag
)
{
	if ( tag && tag->hasOption("search_paths") ) {
		std::string paths = tag->getOption<string>("search_paths");
		utility::vector1< std::string > path_vector = utility::string_split( paths, ',' );
		if ( path_vector.size() != 0 ) {
			search_paths_.clear();
			search_paths_.reserve( path_vector.size() );
			for ( std::string const & path : path_vector ) {
				if ( path.size() == 0 ) {
					search_paths_.push_back( "./" );
				} else if ( path[ path.size() - 1 ] != '/' ) {
					search_paths_.push_back( path + '/' );
				} else {
					search_paths_.push_back( path );
				}
			}
		}
	}
}

void
FileSystemResourceLocator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
)
{
	using namespace utility::tag;
	AttributeList attrs;
	attrs + XMLSchemaAttribute( "search_paths", xs_string, "The comma-separated list of directories where the requested files should be looked for; either an empty string or './' can be used to designate the present working directory" );

	xsd_type_definition_w_attributes( xsd, classname(), "The file system resource locator will interpret the input_id for a Resource as a file name and search"
		" for a file with that name.", attrs );
}


} // namespace locator
} // namespace resource_manager
} // namespace basic
