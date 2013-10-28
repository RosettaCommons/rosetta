// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsExplicitDefiner.cc
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

// Unit Headers
#include <protocols/loops/loops_definers/LoopsExplicitDefiner.hh>
#include <protocols/loops/Loop.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <basic/datacache/DataMap.fwd.hh>


// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>

// Boost Headers
#include <boost/foreach.hpp>

// C++ Headers
#include <string>
#include <utility/excn/Exceptions.hh>
#include <sstream>

/// Macros are not properly caught and passed along by my #inclusion
/// cleanup script
#define foreach BOOST_FOREACH



using std::string;
using std::endl;
using std::stringstream;
using core::Real;
using core::pose::Pose;
using basic::datacache::DataMap;
using utility::tag::TagCOP;
using utility::vector0;
using basic::Tracer;

namespace protocols {
namespace loops {
namespace loops_definers {


static Tracer TR("protocols.loops.loops_definers.LoopsExplicitDefiner");


LoopsExplicitDefiner::LoopsExplicitDefiner() :
	loop_list_()
{}

LoopsExplicitDefiner::~LoopsExplicitDefiner() {}

LoopsExplicitDefiner::LoopsExplicitDefiner(LoopsExplicitDefiner const & src) : LoopsDefiner(src),
		loop_list_(src.loop_list_)
{}

/// @brief Create another loops definer of the type matching the most-derived
/// version of the class.
LoopsDefinerOP
LoopsExplicitDefiner::clone(
) const {
	return new LoopsExplicitDefiner(*this);
}

SerializedLoop
LoopsExplicitDefiner::parse_loop_tag(
	TagCOP const tag,
	string const & loops_name
) {

	SerializedLoop loop;

	if(tag->hasOption("start")){
		loop.start = tag->getOption<Size>("start");
	} else {
		stringstream err_msg;
		err_msg
			<< "Tag " << tag->getName() << " with name "
			<< "'" << loops_name << "' does not have the expected 'start' field." << endl;
		utility_exit_with_message(err_msg.str());
	}

	if(tag->hasOption("stop")){
		loop.stop = tag->getOption<Size>("stop");
	} else {
		stringstream err_msg;
		err_msg
			<< "Tag " << tag->getName() << " with name "
			<< "'" << loops_name << "' does not have the expected 'stop' field." << endl;
		utility_exit_with_message(err_msg.str());
	}

	loop.cut = tag->getOption<Size>("cut", 0);
	loop.skip_rate = tag->getOption<Real>("skip_rate", 0.0);
	loop.extended = tag->getOption<bool>("extended", false);

	return loop;
}


/// @brief Used to parse an xml-like tag to load parameters and properties.
void
LoopsExplicitDefiner::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap const &,
	Pose const &
) {

	if(!tag->hasOption("name")){
		throw utility::excn::EXCN_RosettaScriptsOption(
			"Unable to create unnamed LoopsDefiner (type: Loops)" );
	}
	string const loops_name(tag->getOption<string>("name"));


	vector0< TagCOP >::const_iterator begin=tag->getTags().begin();
	vector0< TagCOP >::const_iterator end=tag->getTags().end();

	for(; begin != end; ++begin){
		TagCOP loop_tag= *begin;

		if(loop_tag->getName() != "loop"){
			TR.Error
				<< "Please include only tags with name 'loop' "
				<< "as subtags of a 'Loops' tag" << endl
				<< "Tag with name '" << loop_tag->getName() << "' is invalid" << endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}

		SerializedLoop loop = parse_loop_tag(loop_tag, loops_name);
		loop_list_.push_back(loop);
	}
}

SerializedLoopList
LoopsExplicitDefiner::apply(
	Pose const &
) {
	return loop_list_;
}


} //namespace
} //namespace
} //namespace
