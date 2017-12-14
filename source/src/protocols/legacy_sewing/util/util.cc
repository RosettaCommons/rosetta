// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <protocols/legacy_sewing/util/util.hh>

//Protocol headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/analysis/LoopAnalyzerMover.hh>
#include <protocols/simple_moves/MinMover.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/legacy_sewing.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>


namespace protocols {
namespace legacy_sewing  {

static basic::Tracer TR( "protocols.legacy_sewing.util" );

std::map<core::id::AtomID, core::id::AtomID>
largest_continuous_atom_map(
	std::map<core::id::AtomID, core::id::AtomID> const & atom_map
) {
	auto it = atom_map.begin();
	auto it_end = atom_map.end();

	//Find the largest continuous stretch of atoms in the alignment
	core::Size prev_resnum_1(0);
	core::Size prev_resnum_2(0);
	std::map<core::id::AtomID, core::id::AtomID> largest_continuous;
	std::map<core::id::AtomID, core::id::AtomID> current_stretch;
	for ( ; it != it_end; ++it ) {
		if ( current_stretch.size() > 0 ) {
			if ( !( (prev_resnum_1 == it->first.rsd() && prev_resnum_2 == it->second.rsd()) ||
					(prev_resnum_1 == it->first.rsd()-1 && prev_resnum_2 == it->second.rsd()-1) )
					) {
				if ( current_stretch.size() > largest_continuous.size() ) {
					largest_continuous = current_stretch;
				}
				current_stretch.clear();
				current_stretch.insert(*it);
			}
		}
		current_stretch.insert(*it);
		prev_resnum_1 = it->first.rsd();
		prev_resnum_2 = it->second.rsd();
	}
	if ( current_stretch.size() > largest_continuous.size() ) {
		largest_continuous = current_stretch;
	}

	return largest_continuous;
}

}
}
