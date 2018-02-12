// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//#include <devel/init.hh>
//#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//#include <utility/exit.hh>
//#include <utility/excn/Exceptions.hh>
//#include <protocols/rosetta_scripts/util.hh>
//#include <fstream>
#include <sstream>

void add_protocols(
	utility::tag::Tag const & protocols_tag,
	std::stringstream & ss
){
	runtime_assert( protocols_tag.getName() == "PROTOCOLS");
	for ( auto const & possible_stage_tag : protocols_tag.getTags() ) {
		if ( possible_stage_tag->getName() != "Stage" ) continue;

		for ( auto const & add_or_sort_tag : possible_stage_tag->getTags() ) {
			if ( add_or_sort_tag->getName() == "Add" ) {
				ss << "\t\t";
				add_or_sort_tag->write( ss , 0 );
				//ss << std::endl;
			} else if ( add_or_sort_tag->getName() == "Sort" ) {
				std::string tag = add_or_sort_tag->to_string();
				auto pos = tag.find("Sort");
				tag.replace( pos, 4, "Add" );
				ss << "\t\t" << tag;// << std::endl;
			} else {
				utility_exit_with_message( "This is neither Add nor Sort: " + add_or_sort_tag->getName() );
			}
		}
	}

}

std::string
reverse_convert(
	utility::tag::Tag const & multistage_rosetta_scripts_tag
){
	runtime_assert( multistage_rosetta_scripts_tag.getName() == "JobDefinitionFile");
	std::stringstream ss;

	ss << "<ROSETTASCRIPTS>" << std::endl;

	for ( auto const & common_or_job_tag : multistage_rosetta_scripts_tag.getTags() ) {
		if ( common_or_job_tag->getName() != "Common" ) continue;

		for ( auto const & subtag : common_or_job_tag->getTags() ) {
			if ( subtag->getName() == "PROTOCOLS" ) {
				ss << "\t<PROTOCOLS>" << std::endl;
				add_protocols( * subtag, ss );
				ss << "\t</PROTOCOLS>" << std::endl;
			} else {
				subtag->write( ss, 1 );
			}
		}
	}

	ss << "</ROSETTASCRIPTS>" << std::endl;
	return ss.str();
}
