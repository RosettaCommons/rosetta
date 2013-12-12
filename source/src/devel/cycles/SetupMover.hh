// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_devel_cycles_SetupMover_hh
#define INCLUDED_devel_cycles_SetupMover_hh

#include <protocols/moves/Mover.hh>

namespace devel {
namespace cycles {

using namespace core;

class SetupMover : public protocols::moves::Mover {
	
	public:
		/// Initialize the mover.
		SetupMover();

		/// Clean up the mover.
		~SetupMover();

		/// Return "SetupMover".
		std::string get_name() const;

		/// Prepare the pose to be used by the other cyclic movers.
		void apply(pose::Pose &pose);

	private:
		/// Remove the free amino and carboxyl groups from the termini.
		void remove_termini_patches(pose::Pose& pose);

		/// Add the amide hydrogen back onto the former N-terminus.
		void restore_amide_hydrogen(pose::Pose& pose);

};

} // End 'cycles' namespace
} // End 'devel' namespace

#endif
