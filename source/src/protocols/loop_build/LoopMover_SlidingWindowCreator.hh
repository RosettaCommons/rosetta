// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loop_build/LoopMover_SlidingWindowCreator.hh
/// @brief  Header for LoopMover_SlidingWindowCreator
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_loop_build_LoopMover_SlidingWindowCreator_hh
#define INCLUDED_protocols_loop_build_LoopMover_SlidingWindowCreator_hh

#include <protocols/moves/MoverCreator.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loop_build {

/// @brief creator for the LoopMover_SlidingWindowCreator class
class LoopMover_SlidingWindowCreator : public moves::MoverCreator {
public:
	LoopMover_SlidingWindowCreator() {}

	virtual ~LoopMover_SlidingWindowCreator();

	moves::MoverOP create_mover() const;

	virtual std::string keyname() const;

};

} //namespace
} //namespace

#endif
