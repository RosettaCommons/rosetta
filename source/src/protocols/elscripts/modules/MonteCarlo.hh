// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/elscripts/modules/MonteCarlo.hh
/// @brief  MonteCarlo module for ElScripts
/// @author Ken Jung

#ifndef INCLUDED_protocols_elscripts_module_MonteCarlo_hh
#define INCLUDED_protocols_elscripts_module_MonteCarlo_hh
#ifdef USELUA
namespace protocols {
namespace elscripts {
namespace modules {

// we're going with the easiest form of lua modules because I don't know what I'm doing
std::string MonteCarlo = 
"	----- MonteCarlo Module ------\n"
"	if modules == nil then\n"
"		modules = {}\n"
"		setmetatable( modules, {__index = _G} )\n"
"	end\n"
"	modules.MonteCarlo = function ( pose, mover, filter, trials, temp)\n"
"		if pose == nil then error(\"Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )\") end\n"
"		if mover == nil then error(\"Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )\") end\n"
"		if filter == nil then error(\"Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )\") end\n"
"		if trials == nil then trials = 10 end\n"
"		if temp == nil then temp = 1 end\n"
"\n"
"		math.randomseed( os.time() )\n"
"		local oldscore = filter( pose )\n"
"		local accepted = pose:clone()\n"
"		local working_pose = pose:clone()\n"
"		local accept_counter = 0\n"
"		for i=1,trials do\n"
"			mover( working_pose )\n"
"			local currentscore = filter( working_pose )\n"
"			if currentscore == nil then error('Filter function returning nil, fix that') end\n"
"			if currentscore < oldscore then \n"
"				accept_counter = accept_counter + 1\n"
"				accepted = working_pose\n"
"				oldscore = currentscore\n"
"			elseif math.random() >= math.exp( math.min(40, math.max( -40, (oldscore - currentscore)/temp ) ) ) then\n"
"				accept_counter = accept_counter + 1\n"
"				accepted = working_pose\n"
"				oldscore = currentscore\n"
"			end\n"
"		end\n"
"		print('Finished Monte Carlo. Out of '..trials..' trials '..accept_counter..' accepted')\n"
"	end\n";

} //modules
} //elscripts
} //protocols
#endif
#endif
