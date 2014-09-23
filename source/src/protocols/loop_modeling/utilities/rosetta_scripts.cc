// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	LoopMoverOP loop_mover = utility::pointer::dynamic_pointer_cast< protocols::loop_modeling::LoopMover > ( base_mover );

	if (loop_mover.get() == NULL) {
		stringstream message;
		message << "<" << tag->getName() << "> is not a loop mover, so ";
		message << "cannot be used in the <LoopModeler> protocol." << endl;
		throw utility::excn::EXCN_Msg_Exception(message.str());
	}

	return loop_mover;
}

}
}
}
