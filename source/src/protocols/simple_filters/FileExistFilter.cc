// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/FileExistFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/FileExistFilter.hh>
#include <protocols/simple_filters/FileExistFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <fstream>
#include <iostream>
namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.FileExistFilter" );

protocols::filters::FilterOP
FileExistFilterCreator::create_filter() const { return protocols::filters::FilterOP( new FileExistFilter ); }

std::string
FileExistFilterCreator::keyname() const { return "FileExist"; }

//default ctor
FileExistFilter::FileExistFilter() :
	protocols::filters::Filter( "FileExist" ),
	filename_( "" ),
	ignore_zero_byte_( false )
{}

FileExistFilter::~FileExistFilter() {}

void
FileExistFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	filename_ = tag->getOption< std::string >( "filename" );
	ignore_zero_byte( tag->getOption< bool >( "ignore_zero_byte", false ) );
}

bool
FileExistFilter::apply( core::pose::Pose const & pose ) const {
	return compute( pose );
}

void
FileExistFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"File "<<filename_<<" exists? " << compute( pose )<<'\n';
}

core::Real
FileExistFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
FileExistFilter::compute(
	core::pose::Pose const & /*pose*/
) const {
	using namespace std;

	ifstream infile;
	infile.open( filename_.c_str(), ios::in );
	if ( !infile.good() ) {
		return false;
	}

	if ( !ignore_zero_byte() ) {
		return true;
	}

	/// if the file is there and we're ignoring zero-byte files, we return true only if the file contains information
	core::Size const begin = infile.tellg();
	infile.seekg( 0, ios::end );
	core::Size const end = infile.tellg();
	return( end - begin > 0 );
}

void FileExistFilter::filename( std::string const & f )
{
	filename_ = f;
}

std::string FileExistFilter::filename() const
{
	return filename_;
}

}
}
