// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Protocol headers
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/moves/MoverFactory.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>

// C++ headers
#include <string>
#include <sstream>

using namespace std;

using utility::tag::TagCOP;
using utility::pointer::dynamic_pointer_cast;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using core::pose::Pose;

namespace protocols {
namespace loop_modeling {
namespace utilities {

LoopMoverOP loop_mover_from_tag(
	TagCOP tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose) {

	using protocols::moves::MoverOP;
	using protocols::moves::MoverFactory;

	MoverOP base_mover = MoverFactory::get_instance()->newMover(
		tag, data, filters, movers, pose);
	LoopMoverOP loop_mover = dynamic_pointer_cast< LoopMover > ( base_mover );

	if ( ! loop_mover ) {
		stringstream message;
		message << "<" << tag->getName() << "> is not a loop mover, so ";
		message << "cannot be used in the <LoopModeler> protocol." << endl;
		throw utility::excn::EXCN_Msg_Exception(message.str());
	}

	return loop_mover;
}

/// @brief return a Loops object from tags
/// @details returns a LoopsOP object from tags "loops_file" or a Loop subtag.  Moved from loop_modeling/LoopMover.cc.
protocols::loops::LoopsOP parse_loops_from_tag( utility::tag::TagCOP tag ) {
	protocols::loops::LoopsOP parsed_loops;

	if ( tag->hasOption("loops_file") ) {
		string loops_file = tag->getOption<string>("loops_file");
		parsed_loops = LoopsOP( new Loops(loops_file) );
	}

	for ( utility::tag::TagCOP subtag : tag->getTags("Loop") ) {
		core::Size start = subtag->getOption<Size>("start");
		core::Size stop = subtag->getOption<Size>("stop");
		core::Size cut = subtag->getOption<Size>("cut", 0);
		core::Real skip_rate = subtag->getOption<Real>("skip_rate", 0.0);
		bool extended = subtag->getOption<bool>("rebuild", false);

		if ( ! parsed_loops ) {
			parsed_loops = LoopsOP( new Loops() );
		}
		parsed_loops->add_loop(start, stop, cut, skip_rate, extended);
	}

	return parsed_loops;  //note may be NULL if no loops were found!
}


}
}
}
