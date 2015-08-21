// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loops_definers/util.cc
/// @brief Utility functions useful in LoopDefiner classes.
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>


#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace loops {
namespace loops_definers {


LoopsOP
load_loop_definitions(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	core::pose::Pose const & pose
) {
	using namespace loops;
	using namespace loops_definers;

	std::string loops_str( tag->getOption< std::string >( "loops" ) );

	if ( data.has("loops_definers", loops_str) ) {
		LoopsDefinerOP loops_definer(
			data.get_ptr< LoopsDefiner >("loops_definers", loops_str));
		return LoopsOP( new Loops(loops_definer->apply(pose)) );
	} else {
		return loops_from_string(loops_str, pose );
	}

	return 0;

}


} // namespace
} // namespace
} // namespace
