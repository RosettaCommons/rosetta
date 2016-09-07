// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief LoopRelaxThreadingMover.hh
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_LoopRelaxThreadingMover_hh
#define INCLUDED_protocols_comparative_modeling_LoopRelaxThreadingMover_hh

#include <core/types.hh>
#include <core/fragment/FragSet.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/comparative_modeling/LoopRelaxThreadingMover.fwd.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

// Really simple class - just takes a job from ThreadingInputter and
// runs looprelax. Needs JD2 !!
class LoopRelaxThreadingMover : public moves::Mover {
public:
	// constructor
	LoopRelaxThreadingMover() {}

	void setup();

	void apply( core::pose::Pose & pose ) override;
	// Undefinded comminting out to fix PyRosetta build  bool apply_mt( core::pose::Pose & pose );
	std::string get_name() const override;

private:
	utility::vector1< core::fragment::FragSetOP > frag_libs_;

	// Read parameters
	core::Size max_loop_rebuild;
	core::Real loop_rebuild_filter;
	std::string remodel;
	std::string relax;
};

} // namespace comparative_modeling
} // namespace protocols

#endif // INCLUDED_protocols_comparative_modeling_LoopRelaxThreadingMover_HH
