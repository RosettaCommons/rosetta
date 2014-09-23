// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/floppy_tail/FloppyTailMover.hh
/// @brief FloppyTail headers
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_floppy_tail_FloppyTailMover_hh
#define INCLUDED_protocols_floppy_tail_FloppyTailMover_hh

// Unit Headers
#include <protocols/floppy_tail/FloppyTailMover.fwd.hh>

// Project Headers

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>

#include <core/fragment/ConstantLengthFragSet.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace floppy_tail {

///@brief FloppyTail mover
class FloppyTailMover : public protocols::moves::Mover {

public:

	///@brief ctor with no arguments
	FloppyTailMover();

	///@brief copy ctor
	FloppyTailMover( FloppyTailMover const & rhs );

	///@brief assignment operator
	FloppyTailMover & operator=( FloppyTailMover const & rhs );

	virtual ~FloppyTailMover();
    
    ///@brief set the movemap instead of initializing it from cmd-line
	virtual void set_movemap(core::kinematics::MoveMapOP const movemap);
    
	virtual void set_fa_scorefxn(core::scoring::ScoreFunctionOP const fa_scorefxn);
    
	virtual void set_cen_scorefxn(core::scoring::ScoreFunctionOP const cen_scorefunction);
    
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const { return "FloppyTailMover"; }

	virtual protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new FloppyTailMover ); }

	virtual bool reinitialize_for_each_job() const { return false; }

	virtual bool reinitialize_for_new_input() const { return true; }

private:
	///@brief init_on_new_input system allows for initializing these details the first time apply() is called.  The job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply()).
	virtual void init_on_new_input(core::pose::Pose const & pose);

private:
	core::Size start_;
	core::Size stop_;
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP centroid_scorefunction_;
	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP movemap_;
	core::kinematics::MoveMapOP movemap_lesstail_;
	///@brief stored so that it can be generated in the init_on_new_input function
	core::kinematics::FoldTreeOP foldtree_;
	core::fragment::ConstantLengthFragSetOP fragset3mer_;

};

} //floppy_tail
} //protocols

#endif //INCLUDED_protocols_floppy_tail_FloppyTailMover_hh
