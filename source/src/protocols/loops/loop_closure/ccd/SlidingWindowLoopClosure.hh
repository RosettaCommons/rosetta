// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/SlidingWindowLoopClosure.hh
/// @brief header file for SlidingWindowLoopClosure protocol
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange
/// @author James Thompson


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_SlidingWindowLoopClosure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_SlidingWindowLoopClosure_hh

// Unit Headers
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/fragment/FragID.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/OrderedFragSet.hh>

//#include <protocols/simple_moves/FragmentMover.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>

#include <protocols/constraints_additional/ConstraintEvaluator.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

extern std::string const VDW_FRAG_STORE;
extern std::string const SCORE_FRAG_STORE;
extern std::string const RMSD_FRAG_STORE;

class SlidingWindowLoopClosure : public moves::Mover {
public:
   /// @brief constructor: supply fragsets for fragment moves
	SlidingWindowLoopClosure(
			core::fragment::FragSetCOP fragset,
			core::scoring::ScoreFunctionOP scorefxn,
			core::kinematics::MoveMapCOP movemap
	);

	//@brief just set defaults -- expects fragset, scorefxn and movemap to be set later
	SlidingWindowLoopClosure();

	~SlidingWindowLoopClosure();

	using moves::Mover::apply;

	//@brief run find fragments that close loop  (if ideal loop closing: such that the less_cut pose is close RMSD <0.1 to pose more_cut)
	// returns less_cut and more_cut with best fragment already applied..
	virtual void apply( core::pose::Pose& more_cut, core::pose::Pose& less_cut );
	virtual std::string get_name() const;

	//@brief run find fragments that close loop  (if ideal loop closing: such that the less_cut pose is close RMSD <0.1 to pose more_cut)
	// returns less_cut and more_cut with best fragment already applied..
	virtual void sample_loops( core::pose::Pose& more_cut, core::pose::Pose& less_cut );

	//@brief run find fragments that close loop  (if ideal loop closing: such that the less_cut pose is close RMSD <0.1 to pose more_cut)
	// returns less_cut and more_cut with best fragment already applied..
	virtual void select_final_loop( core::pose::Pose& more_cut, core::pose::Pose& less_cut );

	// returns less_cut and more_cut with best fragment already applied..
	virtual void apply( core::pose::Pose& more_cut );

	static void setPoseExtraScore( core::pose::Pose &pose );

	//@brief return the list of collected fragments -- here it actually will be more than one Frame
	core::fragment::FragSetCOP
	closure_fragments() const {
		return closure_fragments_;
	}

	void keep_fragments( bool setting = true ) {
		bKeepFragments_ = setting;
	}

	void output_debug_structure( core::pose::Pose & pose, std::string prefix );


	//@brief returns current movemap
	core::kinematics::MoveMap const & movemap() const {
		return *movemap_;
	}

	//@brief sets the movemap
	void movemap( core::kinematics::MoveMapCOP movemap );

	core::scoring::ScoreFunction const & scorefxn() const {
	 return *scorefxn_;
	}

	void scorefxn( core::scoring::ScoreFunctionOP sfxn ) {
	 scorefxn_ = sfxn;
	}

	//@brief set fragments for loop-sampling
	void
	fragments( core::fragment::FragSetCOP frags ) {
		ss_info_ = core::fragment::SecondaryStructureOP( new core::fragment::SecondaryStructure( *frags ) );
		fragset_ = frags;
	}


	void scored_frag_cycle_ratio( core::Real setting ) {
		scored_frag_cycle_ratio_ = setting;
	}

	void short_frag_cycle_ratio ( core::Real setting ) {
		short_frag_cycle_ratio_ = setting;
	}
	//   //@brief override cycle setting
	//   void set_cycles( core::Real cycle_ratio = 1.0 );
	void set_bIdealLoopClosing( bool setting ) {
		bIdealLoopClosing_ = setting;
	}

	bool bIdealLoopClosing() const {
		return bIdealLoopClosing_;
	}

	void set_chainbreak_max( core::Real setting ) {
		chainbreak_max_ = setting;
	}

	void
	set_evaluation( evaluation::MetaPoseEvaluatorOP ev ) { 	evaluator_ = ev; };

	void
	set_loop( Loop const& loop_in ) {
		loop_ = loop_in;
	}

	Loop determine_loop( core::pose::Pose const& more_cut, core::pose::Pose & less_cut );

protected:
  typedef std::list< std::pair< core::Real, Loop > > WindowList;
  void generate_window_list( Size loop_size, WindowList& window_list ) const;


	void set_defaults();

	/// @details also modifies the internal scorefxn_ variable
	/// (But not necessarily to same function as is returned.)
	core::scoring::ScoreFunctionOP setup_frag_scorefxn();

	/// @brief process fragments for loop-closing prowess, keep track of best_fragment
	///  return nr_of new good_loops ( vdw criterion )
	Size
	process_fragments(
    core::fragment::FrameList& frame_in,
		core::pose::Pose const& more_cut,
		core::pose::Pose const& loop_pose
	);

//   /// @brief replace scorefxn
//   void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ) {
//     scorefxn_ = scorefxn;
//   }

//   void set_loop( Loop const& loop_in ) {
//     loop_ = loop_in;
//   }

//   Loop const& loop() const {
//     return loop_;
//   }

//   void init_mc();

//   void set_movemap( core::kinematics::MoveMapCOP mm ) {
//     movemap_ = mm;
//   }

//   void set_fragset( core::fragment::FragSetCOP frags ) {
//     fragset_ = frags;
//   }

	bool bQuickTest() const {
		return bQuickTest_;
	}

	/// return the score used for filtering ( scorefxn_ + filter_cst_weight_ * filter_cst_ )
	core::Real filter_score( core::pose::Pose& pose );


protected:
	// min_begin and max_end  as well as cutpoint are stored here after call to init
	Loop loop_;

	core::scoring::ScoreFunctionOP scorefxn_;  //score3

	//@brief movemap --> which dofs can be moved during loops
	core::kinematics::MoveMapCOP movemap_;

	//@brief a MonteCarlo object -- set_default_mc() , access: mc()
	//moves::MonteCarloOP mc_; not used

	core::fragment::FragSetCOP fragset_;

	core::fragment::SecondaryStructureOP ss_info_;

	core::Size min_loop_size_; // = 6;
	core::Size max_loop_size_; // = 12;
	core::Size min_good_loops_; // = 3;
	core::Size min_breakout_good_loops_; //=5;
	core::Real vdw_delta_; // = 0.5;
	core::Real chainbreak_max_; //=2.0
	core::Real score_delta_; //
	core::Real scored_frag_cycle_ratio_;
	core::Real short_frag_cycle_ratio_;
	core::scoring::ScoreType vdw_score_type_; // = scoring::vdw

	bool bKeepFragments_;
	core::fragment::OrderedFragSetOP closure_fragments_;

	core::Real best_score_;
	core::fragment::FragID best_fragment_;

	bool bQuickTest_;

	//@brief if this flag is set the changes are made to the pose without a chainbreak ( new fold-tree )
	// thus the loop-fragments are applied to either an idealized chain (if true) or to a loop with chainbreak
	// loop-fragments are ranked by score of the resulting pose
	bool bIdealLoopClosing_; //default=true,

	evaluation::MetaPoseEvaluatorOP evaluator_;

	constraints_additional::ConstraintEvaluatorOP filter_cst_;
	core::Real filter_cst_weight_;

};

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_SlidingWindowLoopClosure_hh
