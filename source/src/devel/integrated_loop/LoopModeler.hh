// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman


#ifndef INCLUDED_devel_integrated_loop_LoopModeler_hh
#define INCLUDED_devel_integrated_loop_LoopModeler_hh

//Movers
#include <protocols/moves/Mover.hh>

//devel headers
#include <devel/integrated_loop/LoopManager.hh>

//core headers
#include <core/types.hh>

//Loops
#include <protocols/loops/LoopClass.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/loops/looprelax_protocols.hh>

//Fragments
#include <core/fragment/ConstantLengthFragSet.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

///////////////////////////////////////////
///
///  base class for loop modeling and loop refining
///
///////////////////////////////////////////

namespace protocols {
namespace moves {

//v typedef utility::vector1< protocols::Loop > Loops;

class LoopMover : public Mover {

public:
	//@brief constructor
	LoopMover(){}

	//@brief constructor
	LoopMover(
		Loops LoopList,
		core::scoring::ScoreFunctionOP scorefxn
	) : Mover(),
			LoopList_( LoopList ),
			scorefxn_( scorefxn )
	{
		Mover::type("LoopMover");
		set_default();
	}

	protocols::moves::MonteCarloOP get_mc( core::pose::Pose & pose );

protected:
	void apply( core::pose::Pose & pose ) override = 0;
	virtual void apply_mod( core::pose::Pose & pose  ) = 0;
	void set_default();
	void set_default_mc( core::pose::Pose & pose );

	void set_movemap( Loops const & LoopsToPerturb,
		core::kinematics::MoveMapOP movemap_  );
	void set_movemap( protocols::Loop const & ThisLoop,
		core::kinematics::MoveMapOP one_loop_movemap );

	void set_one_loop_fold_tree( core::pose::Pose & pose,
		protocols::Loop const & ThisLoop );
	void set_loops_fold_tree( core::pose::Pose & pose,
		Loops const & LoopsToPerturb );

	//vprotected:
	Loops LoopList_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< utility::vector1< Size > > CompositeLoopIndex_;
	core::kinematics::MoveMapOP movemap_;
	Real mc_temperature_;//default temperature
	Size nmoves_;
	protocols::moves::MonteCarloOP mc_;
	protocols::minimization_packing::MinMoverOP min_mover_;
	bool mc_created;

};

/////////////////////////////////////////////////////
///
///  derived class : loop modeling by fragment insertion
///
///////////////////////////////////////////////////////

class LoopModeler : public LoopMover {

public:
	LoopModeler(
		Loops LoopList,
		core::scoring::ScoreFunctionOP scorefxn,
		core::fragment::FragSetCOP fragset3mer,
		core::fragment::FragSetCOP fragset9mer
	) : LoopMover( LoopList, scorefxn),
			fragset3mer_( fragset3mer ),
			fragset9mer_( fragset9mer )
	{
		Mover::type("LoopModeler");
	}

	void apply( core::pose::Pose & pose ) override;
	void apply_mod( core::pose::Pose & pose ) override;

private:
	core::fragment::FragSetCOP fragset3mer_;
	core::fragment::FragSetCOP fragset9mer_;


};


/////////////////////////////////////////////////////
///
///  derived class : loop refinement (no fragment insertions)
///
///////////////////////////////////////////////////////

class LoopRefiner : public LoopMover {

public:
	LoopRefiner(
		Loops LoopList,
		core::scoring::ScoreFunctionOP scorefxn
	) : LoopMover( LoopList, scorefxn )
	{
		Mover::type("LoopRefiner");
	}

	void apply( core::pose::Pose & pose ) override;
	void apply_mod( core::pose::Pose & pose ) override;
private:

	protocols::moves::SequenceMoverOP small_move_rot_trial_mover( core::kinematics::MoveMapOP movemap_one_loop );
	protocols::moves::SequenceMoverOP shear_move_rot_trial_mover( core::kinematics::MoveMapOP movemap_one_loop );
	protocols::minimization_packing::RotamerTrialsMoverOP rotamer_trial_mover( core::pose::Pose & pose );

	//	moves::SequenceMoverOP SmallMovesRotTrial;
	//	moves::SequenceMoverOP ShearMovesRotTrial;
	//	protocols::minimization_packing::RotamerTrialsMoverOP RotTrial;

};


} //moves
}// protocols


#endif
