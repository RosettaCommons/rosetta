// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Creator class for DisulfideInsertionMover
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Mar. 30, 2015

//
// DisulfideInsertionMoverCreator
//

#include <protocols/simple_moves/DisulfideInsertionMoverCreator.hh>
#include <protocols/simple_moves/DisulfideInsertionMover.hh>

namespace protocols {
namespace simple_moves {

protocols::moves::MoverOP
DisulfideInsertionMoverCreator::create_mover() const {
	core::Size const DEFAULT_CHAIN = 1;
	protocols::moves::MoverOP mover( new protocols::simple_moves::DisulfideInsertionMover(DEFAULT_CHAIN) );
	return mover;
}

std::string
DisulfideInsertionMoverCreator::keyname() const {
	return "DisulfideInsertion";
}

} //namespace simple_moves
} //namespace protocols
