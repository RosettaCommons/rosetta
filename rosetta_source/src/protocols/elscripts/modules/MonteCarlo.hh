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

// we're going with the easiest form of lua modueles because I don't know what I'm doing
std::string MonteCarlo = R"DELIM(
	----- MonteCarlo Module ------
	if modules == nil then
		modules = {}
		setmetatable( modules, {__index = _G} )
	end
	modules.MonteCarlo = function ( pose, mover, filter, trials, temp)
		if pose == nil then error("Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )") end
		if mover == nil then error("Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )") end
		if filter == nil then error("Syntax: modules.MonteCarlo( pose, mover, filter, optional: trials, temp )") end
		if trials == nil then trials = 10 end
		if temp == nil then temp = 1 end

		math.randomseed( os.time() )
		local oldscore = filter( pose )
		local accepted = pose:clone()
		local working_pose = pose:clone()
		local accept_counter = 0
		for i=1,trials do
			mover( working_pose )
			local currentscore = filter( working_pose )
			if currentscore == nil then error('Filter function returning nil, fix that') end
			if currentscore < oldscore then 
				accept_counter = accept_counter + 1
				accepted = working_pose
				oldscore = currentscore
			elseif math.random() >= math.exp( math.min(40, math.max( -40, (oldscore - currentscore)/temp ) ) ) then
				accept_counter = accept_counter + 1
				accepted = working_pose
				oldscore = currentscore
			end
		end
		print('Finished Monte Carlo. Out of '..trials..' trials '..accept_counter..' accepted')
	end
	)DELIM";

} //modules
} //elscripts
} //protocols
#endif
#endif
