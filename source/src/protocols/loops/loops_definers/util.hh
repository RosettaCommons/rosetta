// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/loops_definers/util.hh
/// @brief Utility functions useful in LoopDefiner classes.
/// @authors Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_util_HH
#define INCLUDED_protocols_loops_loops_definers_util_HH

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loops {
namespace loops_definers {

LoopsOP
load_loop_definitions(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap const & data,
	core::pose::Pose const & pose);


} // namespace
} // namespace
} // namespace

#endif  // include guard
