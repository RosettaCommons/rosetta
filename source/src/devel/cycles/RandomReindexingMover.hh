// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_devel_cycles_RandomReindexingMover_hh
#define INCLUDED_devel_cycles_RandomReindexingMover_hh

#include <protocols/moves/Mover.hh>

namespace devel {
namespace cycles {

using namespace core;

class RandomReindexingMover : public protocols::moves::Mover {

	public:
		/// Initialize the mover.
		RandomReindexingMover() {};

		/// Clean up the mover.
		~RandomReindexingMover() {};

		/// Cyclically reindex the given pose by a random amount.
		void apply(pose::Pose &pose);

		/// Return "RandomReindexingMover".
		std::string get_name() const { return "RandomReindexingMover"; }
};

} // End 'cycles' namespace
} // End 'devel' namespace

#endif
