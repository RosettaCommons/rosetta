// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_rbsegment_relax_AutoRBRelaxMover_hh
#define INCLUDED_protocols_rbsegment_relax_AutoRBRelaxMover_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/rbsegment_relax/RBSegment.fwd.hh>
#include <utility/vector1.hh>

#ifdef PYROSETTA
#include <protocols/rbsegment_relax/RBSegment.hh>
#endif


namespace protocols {
namespace rbsegment_relax {

class AutoRBMover : public protocols::moves::Mover {
public:
	AutoRBMover();


	void apply( core::pose::Pose & pose ) override;


	std::string get_name() const override { return "AutoRBMover"; }

	virtual
	void setup_topology( core::pose::Pose & pose );

private:
	void grow_flexible( core::Size maxlen, core::Size nres, core::Size minlen=0 );

	core::Size nouter_cycles_, ninner_cycles_;

	core::Size loop_melt_;

	core::pack::task::TaskFactoryOP tf_;
	core::scoring::ScoreFunctionOP scorefxn_, fa_scorefxn_;
	core::kinematics::MoveMapOP movemap_;
	utility::vector1< RBSegment > rigid_segs_;
	utility::vector1< RBSegment > rb_chunks_;
	utility::vector1< core::Size > jumps_;
	protocols::loops::Loops loops_;

	std::vector< core::fragment::FragSetOP > frag_libs_;

	bool allowSeqShiftMoves_, allowSSFragInserts_;
};

}
}

#endif
