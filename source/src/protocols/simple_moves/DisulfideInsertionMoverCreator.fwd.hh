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


#ifndef INCLUDED_protocols_simple_moves_DisulfideInsertionMoverCreator_fwd_hh
#define INCLUDED_protocols_simple_moves_DisulfideInsertionMoverCreator_fwd_hh

// Project headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace simple_moves {

class DisulfideInsertionMoverCreator;

typedef utility::pointer::shared_ptr< DisulfideInsertionMoverCreator >  DisulfideInsertionMoverCreatorOP;
typedef utility::pointer::shared_ptr< DisulfideInsertionMoverCreator const >  DisulfideInsertionMoverCreatorCOP;

} // namespace simple_moves
} // namespace protocols

#endif
// INCLUDED_protocols_simple_moves_DisulfideInsertionMoverCreator_fwd_hh
