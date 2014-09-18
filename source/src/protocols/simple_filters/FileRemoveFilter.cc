// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <boost/foreach.hpp>
#include <stdio.h>
#include <fstream>


namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.simple_filters.FileRemoveFilter" );

protocols::filters::FilterOP
FileRemoveFilterCreator::create_filter() const { return new FileRemoveFilter; }

std::string
FileRemoveFilterCreator::keyname() const { return "FileRemove"; }

//default ctor
FileRemoveFilter::FileRemoveFilter() :
protocols::filters::Filter( "FileRemove" ),
delete_content_only_( false )
{
	file_names_.clear();
}

FileRemoveFilter::~FileRemoveFilter() {}

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
	BOOST_FOREACH( std::string const f, file_names_ ){
		if( remove( f.c_str() ) )
			TR<<"Successfully removed "<<f<<std::endl;
		else
			TR<<"File "<<f<<" not found."<<std::endl;
		if( delete_content_only() ){
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
FileRemoveFilter::file_names( utility::vector1< std::string > const f ){
	file_names_ = f;
}

}
}
