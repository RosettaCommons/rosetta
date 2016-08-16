// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author J. Karanicolas

#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_Backrub_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_Backrub_hh

#include <protocols/loops/loop_mover/refine/LoopMover_Backrub.fwd.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
//#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class LoopMover_Refine_Backrub: public LoopMover {
public:
	LoopMover_Refine_Backrub();

	LoopMover_Refine_Backrub(
		protocols::loops::LoopsOP  loops_in
	);

	LoopMover_Refine_Backrub(
		protocols::loops::LoopsOP  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	//destructor
	~LoopMover_Refine_Backrub();

	virtual std::string get_name() const;

	void set_default_settings(){
		redesign_loop = false;
	}

	void set_redesign_loop( bool value = true ){ redesign_loop = value; }
	bool get_redesign_loop(){ return redesign_loop; }

	void set_task_factory( core::pack::task::TaskFactoryOP value );
	bool get_task_factory();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return protocols::moves::MoverOP( new LoopMover_Refine_Backrub(*this) ); // <--- TaskFactory.hh has to be #included here because this class's copy constructor is undefined, and this function is being invoked in the header
	}

	void apply( core::pose::Pose & pose );

protected:

	core::pack::task::TaskFactoryOP task_factory;
	bool redesign_loop;
	virtual basic::Tracer & tr() const;
};

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_Backrub_hh
