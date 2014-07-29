// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Test headers
#include <cxxtest/TestSuite.h>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>
#include <sstream>

template <class MoverSubclass>
utility::pointer::owning_ptr<MoverSubclass> parse_tag(std::string tag_string) {
	std::istringstream tag_stream(tag_string);
	utility::tag::TagCOP tag = utility::tag::Tag::create(tag_stream);
	basic::datacache::DataMap data;
	protocols::filters::Filters_map filters;
	protocols::moves::Movers_map movers;
	core::pose::Pose pose;

	protocols::moves::MoverOP base_mover =
		protocols::moves::MoverFactory::get_instance()->newMover(
				tag, data, filters, movers, pose);
	utility::pointer::owning_ptr<MoverSubclass> mover = 
		dynamic_cast<MoverSubclass*>(base_mover.get());

	TSM_ASSERT("Instantiated the wrong type of mover", mover.get());

	return mover;
}
