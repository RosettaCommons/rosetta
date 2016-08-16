// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_build_LoopBuildMover_hh
#define INCLUDED_protocols_loop_build_LoopBuildMover_hh

#include <core/pose/Pose.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <protocols/loop_build/LoopBuildMover.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace loop_build {


class LoopBuildMover : public moves::Mover {
public:

	LoopBuildMover(protocols::comparative_modeling::LoopRelaxMover loop_relax_mover);

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:
	protocols::comparative_modeling::LoopRelaxMover loop_relax_mover_;

private:
	void setup_loop_definition();
};


} // namespace loop_build
} // namespace protocols

#endif /* INCLUDED_protocols_loop_build_LoopBuildMover_hh */
