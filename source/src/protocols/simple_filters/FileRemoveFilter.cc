// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/FileRemoveFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/FileRemoveFilter.hh>
#include <protocols/simple_filters/FileRemoveFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <cstdio>
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.FileRemoveFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP FileRemoveFilterCreator::create_filter() const { return protocols::filters::FilterOP( new FileRemoveFilter ); }

// XRW TEMP std::string
// XRW TEMP FileRemoveFilterCreator::keyname() const { return "FileRemove"; }

//default ctor
FileRemoveFilter::FileRemoveFilter() :
	protocols::filters::Filter( "FileRemove" ),
	delete_content_only_( false )
{
	file_names_.clear();
}

FileRemoveFilter::~FileRemoveFilter() = default;

void
FileRemoveFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	std::string s;
	s = tag->getOption< std::string >( "filenames" );
	file_names( utility::string_split( s, ',' ) );
	delete_content_only( tag->getOption< bool >( "delete_content_only", false ) );
}

bool
FileRemoveFilter::apply( core::pose::Pose const & ) const {
	using namespace std;
	for ( std::string const & f : file_names_ ) {
		if ( remove( f.c_str() ) ) {
			TR<<"Successfully removed "<<f<<std::endl;
		} else {
			TR<<"File "<<f<<" not found."<<std::endl;
		}
		if ( delete_content_only() ) {
			TR<<"Leaving 0b placeholder for file "<<f<<std::endl;
			ofstream outfile;
			outfile.open( f.c_str(), ios::trunc );
			outfile.close();
		}
	}
	return true;
}

void
FileRemoveFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<compute( pose )<<std::endl;
}

core::Real
FileRemoveFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
FileRemoveFilter::compute(
	core::pose::Pose const &
) const {
	return 1;
}

utility::vector1< std::string >
FileRemoveFilter::file_names() const{
	return file_names_;
}

void
FileRemoveFilter::file_names( utility::vector1< std::string > const & f ){
	file_names_ = f;
}

std::string FileRemoveFilter::name() const {
	return class_name();
}

std::string FileRemoveFilter::class_name() {
	return "FileRemove";
}

void FileRemoveFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("filenames", xs_string, "list of file names separated by comma, e.g., 3r2x_0001.pdb,3r2x_0002.pdb")
		+ XMLSchemaAttribute::attribute_w_default("delete_content_only", xsct_rosetta_bool, "if true, only eliminates the contents of the file but leaves a placeholder file of size 0bytes.", "false");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Remove a file from disk.", attlist );
}

std::string FileRemoveFilterCreator::keyname() const {
	return FileRemoveFilter::class_name();
}

protocols::filters::FilterOP
FileRemoveFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new FileRemoveFilter );
}

void FileRemoveFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FileRemoveFilter::provide_xml_schema( xsd );
}


}
}
