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

#ifndef INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_hh
#define INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_hh

#include <protocols/loops/loop_mover/refine/LoopMover_CCD.fwd.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>

#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

class LoopMover_Refine_CCD: public loop_mover::LoopMover {
public:

	// empty constructor
	LoopMover_Refine_CCD();

	//constructor
	LoopMover_Refine_CCD(
		protocols::loops::LoopsOP  loops_in
	);

	LoopMover_Refine_CCD(
		protocols::loops::LoopsOP  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	//destructor
	~LoopMover_Refine_CCD();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	void set_default_settings();

	void set_redesign_loop( bool value = true ){ redesign_loop_ = value; }
	bool get_redesign_loop(){ return redesign_loop_; }

	void set_task_factory( core::pack::task::TaskFactoryCOP task_factory_in );

	core::pack::task::TaskFactoryCOP get_task_factory() const;

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void outer_cycles( core::Size value ) { outer_cycles_ = value; }
	void max_inner_cycles( core::Size value ) { max_inner_cycles_ = value; }
	void repack_period( core::Size value ) { repack_period_ = value; }
	void temp_initial( core::Real value ) { temp_initial_ = value; }
	void temp_final( core::Real value ) { temp_final_ = value; }
	void set_fold_tree_from_loops( bool const s ){ set_fold_tree_from_loops_ = s; }
	bool set_fold_tree_from_loops() const{ return set_fold_tree_from_loops_; }
	core::kinematics::MoveMapOP move_map() const;
	void move_map( core::kinematics::MoveMapOP mm );
	void set_flank_residue_min(bool value) {flank_residue_min_ = value;} // by JQX
	friend std::ostream &operator<< ( std::ostream &os, LoopMover_Refine_CCD const &mover )
	{
		moves::operator<<(os, mover);
		os << "Outer cycles: " << mover.outer_cycles_ << ",  Max inner cycles: " << mover.max_inner_cycles_ << ",  Repack period: " << mover.repack_period_ << std::endl <<
				"Initial temperature: " << mover.temp_initial_ << ",  Final temperature: " << mover.temp_final_ << std::endl <<
				"Set fold tree from loop?: " << mover.set_fold_tree_from_loops_ << std::endl;
		return os;
	}

protected:
	void read_options();

	/// @brief setup an appropriate movemap for the given loops
	/// @param[in] loops The loops to model.
	/// @param[in] allow_repack Indicates whether or not to allow a position to
	///  repack.
	/// @param[out] movemap Output movemap, all settings added here.
	/// @remarks will enforce the false movemap
	void setup_movemap(
		core::pose::Pose const & pose,
		protocols::loops::Loops const & loops,
		utility::vector1< bool > const & allow_repack,
		core::kinematics::MoveMapOP & movemap
	);

	core::pack::task::TaskFactoryOP task_factory_;
	bool redesign_loop_;
    virtual basic::Tracer & tr() const;

private:
	// parameters with local defaults
	// Overriden by options if specified by user (do not use option defaults), or via setter methods
	core::Size outer_cycles_, max_inner_cycles_, repack_period_;
	core::Real temp_initial_, temp_final_;
	bool packing_isolated_to_active_loops_;
	bool set_fold_tree_from_loops_;
	core::kinematics::MoveMapOP move_map_;
	bool flank_residue_min_; //JQX

}; // LoopMover_Refine_CCD

} //namespace refine
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_refine_LoopMover_CCD_HH
