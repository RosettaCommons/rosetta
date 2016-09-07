// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/abinitio/AbrelaxMover.cc
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)
/// @brief This class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()

#ifndef INCLUDED_protocols_abinitio_AbrelaxMover_hh
#define INCLUDED_protocols_abinitio_AbrelaxMover_hh

// Unit Headers
#include <protocols/abinitio/AbrelaxMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.fwd.hh>
#include <protocols/relax/RelaxProtocolBase.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyBroker.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// C++ headers
#include <string>

#include <protocols/abinitio/FragmentSampler.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class AbrelaxMover : public moves::Mover {
public:
	AbrelaxMover();

	~AbrelaxMover() override;

	void apply( core::pose::Pose &pose ) override;

	std::string get_name() const override;

	void set_defaults();

	void clear(); //set's all pointers to zero

	moves::MoverOP fresh_instance() const override { return moves::MoverOP( new AbrelaxMover() ); }

	FragmentSamplerOP sampling_protocol();

	relax::RelaxProtocolBaseCOP relax_protocol() const;

	relax::RelaxProtocolBaseOP relax_protocol();

	loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol();

	void sampling_protocol( FragmentSamplerOP set);

	void relax_protocol( relax::RelaxProtocolBaseOP set );

	void closure_protocol( loops::loop_closure::ccd::SlidingWindowLoopClosureOP set );

	void post_loop_closure_protocol( moves::MoverOP move );  //e.g. to idealize after loop-closing

	void pre_loop_closure_protocol( moves::MoverOP move ); //e.g. to idealize after loop-closing

	topology_broker::TopologyBrokerOP topology_broker();

private:
	void close_with_idealization( core::pose::Pose &pose);

	topology_broker::TopologyBrokerOP topology_broker_;

	FragmentSamplerOP sampling_protocol_;
	loops::loop_closure::ccd::SlidingWindowLoopClosureOP loop_closure_protocol_;
	relax::RelaxProtocolBaseOP relax_protocol_;
	moves::MoverOP post_loop_closure_protocol_;
	moves::MoverOP pre_loop_closure_protocol_;

	bool b_return_unrelaxed_fullatom_;
};

}
}

#endif  // INCLUDED_protocols_abinitio_AbrelaxMover_hh
