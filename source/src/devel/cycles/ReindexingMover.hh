// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_devel_cycles_ReindexingMover_hh
#define INCLUDED_devel_cycles_ReindexingMover_hh

#include <protocols/moves/Mover.hh>

namespace devel {
namespace cycles {

using namespace core;

class ReindexingMover : public protocols::moves::Mover {
	
	public:
		/// Initialize the mover.
		ReindexingMover(Size offset=0);

		/// Clean up the mover.
		~ReindexingMover() override;

		/// Cyclically reindex all of the residues in the given pose.
		void apply(pose::Pose &pose) override;

		/// Return "ReindexingMover".
		std::string get_name() const override { return "ReindexingMover"; }

		/// Specify how many residues to rotate on future apply() calls.
		void set_offset(Size offset) { offset_ = offset; }

	private:
		int offset_;
};

} // End 'cycles' namespace
} // End 'devel' namespace

#endif
