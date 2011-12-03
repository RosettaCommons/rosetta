// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Mike Tyka
/// @author Daniel J. Mandell

#ifndef INCLUDED_protocols_loops_LoopMover_KIC_hh
#define INCLUDED_protocols_loops_LoopMover_KIC_hh


#include <protocols/loops/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <basic/Tracer.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// C++ Headers


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {

class LoopMover_Perturb_KIC: public IndependentLoopMover {
public:

	LoopMover_Perturb_KIC();

	LoopMover_Perturb_KIC(
		protocols::loops::Loops  loops_in
	);

	LoopMover_Perturb_KIC(
		protocols::loops::Loops  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	//destructor
	~LoopMover_Perturb_KIC();

	virtual std::string get_name() const;

	virtual void  set_extended_torsions(
		core::pose::Pose & pose,
		Loop const & loop
	);

	void set_default_settings();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	LoopResult model_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop
	);

private:

	core::Size max_seglen_; // maximum KIC segment length
	bool recover_low_;
	core::Size max_kic_build_attempts_;
	core::Size remodel_kic_attempts_;
};


class LoopMover_Refine_KIC: public LoopMover {
public:
	//constructors
	LoopMover_Refine_KIC(
		protocols::loops::Loops  loops_in
	);

	LoopMover_Refine_KIC(
		protocols::loops::Loops const loops_in,
		core::scoring::ScoreFunctionCOP  scorefxn
	);

	//destructor
	~LoopMover_Refine_KIC();

	void set_default_settings();

	void set_redesign_loop( bool value = true ){ redesign_loop_ = value; }
	bool get_redesign_loop(){ return redesign_loop_; }

	void set_task_factory( core::pack::task::TaskFactoryOP value );
	bool get_task_factory();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	///@brief update the vector of movemaps, one for each loop in loops_
	void
	update_movemap_vectors(
		core::pose::Pose & pose,
		utility::vector1<core::kinematics::MoveMap> & move_maps );

	///@brief update the vector of vectors of moveable side-chain positions, one for each loop in loops_
	void
	update_allow_sc_vectors(
		core::pose::Pose & pose,
		utility::vector1< utility::vector1< bool > > & allow_sc_vectors );

	void
	set_rottrials_from_kic_segment(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP & rottrials_packer_task,
		Size kic_start,
		Size kic_end );

	void
	set_movemap_from_kic_segment(
		core::pose::Pose & pose,
		core::kinematics::MoveMap & cur_mm,
		Size kic_start,
		Size kic_end );

protected:

	core::pack::task::TaskFactoryOP task_factory;
	bool redesign_loop_;

private:

	core::Real neighbor_dist_; // CB distance to loop to consider scaffold side-chain for rot trials, repack, dfpmin
	core::Size max_seglen_; // maximum KIC segment length
	bool recover_low_; // flag to recover the lowest energy MC conformation rather than the last conformation
	bool min_after_repack_; // should inner cycle repacking steps be followed by minimization
	bool fix_natsc_; // should side-chains neighboring the loop be fixed
	bool optimize_only_kic_region_sidechains_after_move_; // Should we perform rotamer trials and minimization after every
													      // KIC move but only within the neighbor_dist of the KIC segment
};



} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopMover_KIC_HH
