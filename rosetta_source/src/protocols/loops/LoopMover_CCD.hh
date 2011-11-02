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

#ifndef INCLUDED_protocols_loops_LoopMover_CCD_hh
#define INCLUDED_protocols_loops_LoopMover_CCD_hh


#include <protocols/loops/IndependentLoopMover.hh>
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


// C++ Headers

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {

////////////////////////////////////////////////////////////////////////////////////////
/// @details the main function to model one single loop in centroid mode. The
/// modeling algorithm is fragment_ccd_min_trial, which consists of perturbing
/// the loop conformation by fragment insertion , then close the loop by CCD loop
/// closure, then minimize the loop conformation and finally subject it to Monte
/// Carlo acceptance or rejection. The loop conformation will be initialized as
/// extended conformation if it is specified in the loop definition, resembling
/// ab initio loop modeling in real practice. The loop has to be long enough for
/// inserting certain length of fragments.
/////////////////////////////////////////////////////////////////////////////////////////
class LoopMover_Perturb_CCD: public IndependentLoopMover {
public:
	//constructor
	LoopMover_Perturb_CCD(
		protocols::loops::Loops  loops_in
	);

	LoopMover_Perturb_CCD(
		protocols::loops::Loops  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	LoopMover_Perturb_CCD(
		protocols::loops::Loops  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn,
		core::fragment::FragSetOP fragset
	);

	//destructor
	~LoopMover_Perturb_CCD();

	virtual std::string get_name() const;

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	void set_default_settings(){}

	LoopResult model_loop( core::pose::Pose & pose,
	                 protocols::loops::Loop const & loop );

protected:
	std::vector< core::fragment::FragSetOP > frag_libs_;

};

class LoopMover_Refine_CCD: public LoopMover {
public:

	// empty constructor
	LoopMover_Refine_CCD();

	//constructor
	LoopMover_Refine_CCD(
		protocols::loops::Loops  loops_in
	);

	LoopMover_Refine_CCD(
		protocols::loops::Loops  loops_in,
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
		core::kinematics::MoveMap & movemap
	);

	core::pack::task::TaskFactoryOP task_factory_;
	bool redesign_loop_;

private:
	// parameters with local defaults
	// Overriden by options if specified by user (do not use option defaults), or via setter methods
	core::Size outer_cycles_, max_inner_cycles_, repack_period_;
	core::Real temp_initial_, temp_final_;
	bool packing_isolated_to_active_loops_;
	bool set_fold_tree_from_loops_;

}; // LoopMover_Refine_CCD

} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopMover_CCD_HH
