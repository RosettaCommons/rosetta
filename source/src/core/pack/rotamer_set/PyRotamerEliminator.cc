// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/PyRotamerEliminator.cc
/// @brief  Wrappers that allow the user to eliminate rotamers from pyrosetta
/// @author Jack Maguire, jackmaguire1444@gmail.com


// Unit Headers
#include <core/pack/rotamer_set/PyRotamerEliminator.hh>

#include <utility/vector1.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>

namespace core {
namespace pack {
namespace rotamer_set {

RotamerSetsOperationOP
PyRotamerEliminator::clone() const {
	return utility::pointer::make_shared< PyRotamerEliminator >( *this );
}

void
PyRotamerEliminator::alter_rotamer_sets(
	pose::Pose const & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTask const & ptask,
	utility::graph::GraphCOP packer_neighbor_graph,
	RotamerSets & rotamer_sets
) {
	// Performance Warning!
	// Python currently creates copies of everything that is not wrapped in an OP or COP
	// Not only is this slow, but any changes made to rotamer_sets will be made to the copy,
	//     hence the vector-of-bool return value instead of editing rotamer_sets directly
	auto const & rotamers_to_delete =
		function_( pose, sfxn.clone(), ptask.clone(), packer_neighbor_graph, rotamer_sets );

	runtime_assert( rotamers_to_delete.size() == rotamer_sets.nrotamers() );

	rotamer_sets.drop_rotamers( rotamers_to_delete );
}



}
}
}

