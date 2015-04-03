// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/matdes/MatDesGreedyOptMutationMoverMatDesGreedyOptMutationMoverCreator.hh
/// @brief	This is a modified version of Chris King's GreedyOptMutationMover with additional functionality that is currently not compatible with all of the ParetoOpt functionality. Please note that this has been checked into master in its current state in response to requests from others to use this modified version of Chris King's GreedyOptMutationMover. Although this is still a somewhat developmental piece of code, it has currently been left in src/protocols/matdes/ to avoid issues with intra-library level dependencies.  
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_MatDesGreedyOptMutationMoverCreator_hh
#define INCLUDED_protocols_matdes_MatDesGreedyOptMutationMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace matdes {

class MatDesGreedyOptMutationMoverCreator : public protocols::moves::MoverCreator {
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
};

}
}

#endif

