// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/floppy_tail/TailSegmentMover.hh
/// @brief FloppyTail headers
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_devel_splice_TailSegmentMover_hh
#define INCLUDED_devel_splice_TailSegmentMover_hh

// Unit Headers
#include <devel/splice/TailSegmentMover.fwd.hh>

// Project Headers

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>

#include <core/fragment/ConstantLengthFragSet.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>

namespace devel {
namespace splice {

/// @brief FloppyTail mover

class TailSegmentMover : public protocols::moves::Mover {

public:
	/// @brief ctor with no arguments
	TailSegmentMover();

	/// @brief copy ctor
	TailSegmentMover( TailSegmentMover const & rhs );

	/// @brief assignment operator
	TailSegmentMover & operator=( TailSegmentMover const & rhs );

	virtual ~TailSegmentMover();

	/// @brief set the movemap instead of initializing it from cmd-line
	virtual void set_movemap(core::kinematics::MoveMapOP const movemap);

	virtual void set_fa_scorefxn(core::scoring::ScoreFunctionOP const fa_scorefxn);


	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const { return "TailSegmentMover"; }

	virtual protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new TailSegmentMover ); }

	virtual bool reinitialize_for_each_job() const { return false; }

	virtual bool reinitialize_for_new_input() const { return true; }

	void temp_initial( core::Real value ) { temp_initial_ = value; }
	void temp_final( core::Real value ) { temp_final_ = value; }

	void set_task_factory( core::pack::task::TaskFactoryCOP task_factory_in );


private:
	core::Size start_;
	core::Size stop_;
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP centroid_scorefunction_;
	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP movemap_;
	core::kinematics::MoveMapOP movemap_lesstail_;
	/// @brief stored so that it can be generated in the init_on_new_input function
	core::kinematics::FoldTreeOP foldtree_;
	core::Real temp_initial_, temp_final_;


};

} //floppy_tail
} //protocols

#endif //INCLUDED_protocols_floppy_tail_FloppyTailMover_hh
