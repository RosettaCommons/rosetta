// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/elscripts/util.hh
/// @brief Utility functions useful in elscripts
/// @author ken Jung

#ifndef INCLUDED_protocols_elscripts_util_hh
#define INCLUDED_protocols_elscripts_util_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/lua/LuaObject.hh>
#include <utility/lua/LuaIterator.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

namespace protocols {
namespace elscripts {

core::scoring::ScoreFunctionOP parse_scoredef( utility::lua::LuaObject const & scoredef,
	utility::lua::LuaObject const & score_fxns );

void parse_movemapdef( utility::lua::LuaObject const & movemapdef, core::kinematics::MoveMapOP mm );

core::pack::task::TaskFactoryOP
parse_taskdef( utility::lua::LuaObject const & taskdef,
	utility::lua::LuaObject const & tasks );

core::pack::task::TaskFactorySP
sp_parse_taskdef( utility::lua::LuaObject const & taskdef,
	utility::lua::LuaObject const & tasks );

} // protocols
} // elscripts
#endif
