// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsFileDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsFileDefiner.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/Loop.hh>

// Project Headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>


using std::string;
using std::endl;
using std::stringstream;
using utility::tag::TagCOP;
using basic::datacache::DataMap;
using core::pose::Pose;


namespace protocols {
namespace loops {
namespace loops_definers {

LoopsFileDefiner::LoopsFileDefiner() :
	loop_list_()
{}

LoopsFileDefiner::~LoopsFileDefiner() {}

LoopsFileDefiner::LoopsFileDefiner(LoopsFileDefiner const & src) : LoopsDefiner(src),
	loop_list_(src.loop_list_)
{}

/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsFileDefiner::clone(
) const {
	return new LoopsFileDefiner(*this);
}

/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsFileDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &,
	Pose const & p
) {

	if(!tag->hasOption("name")){
		throw utility::excn::EXCN_RosettaScriptsOption(
			"Unable to create unnamed LoopsDefiner (type: LoopsFile)" );
	}
	string const loops_name(tag->getOption<string>("name"));


	string filename;
	if(tag->hasOption("filename")){
		filename = tag->getOption<string>("filename");
	} else {
		stringstream err_msg;
		err_msg << "Tag with name '" << loops_name << "' does not have the expected 'filename' field." << endl;
		throw utility::excn::EXCN_RosettaScriptsOption(err_msg.str());
	}

	LoopsFileIO loops_file_io;
	LoopsFileDataOP lfd = loops_file_io.read_loop_file( filename );
	loop_list_ = lfd->resolve_as_serialized_loops( p );
}

SerializedLoopList
LoopsFileDefiner::apply(
	Pose const &
) {
	return loop_list_;
}

} //namespace
} //namespace
} //namespace
