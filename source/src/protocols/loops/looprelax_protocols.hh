// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopLoopRebuild_protocols.hh
/// @brief
/// @details
///
/// @author James Thompson
/// @author Srivatsan Raman


#ifndef INCLUDED_protocols_loops_looprelax_protocols_hh
#define INCLUDED_protocols_loops_looprelax_protocols_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
//#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <core/fragment/FragSet.fwd.hh>

#ifdef WIN32
#include <core/fragment/FragSet.hh> // WIN32 INCLUDE
#endif

//// C++ headers
#include <string>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {


/// @brief The loop-rebuild protocol
class LoopRebuild: public protocols::moves::Mover {
public:
	/// @brief Construct te protocol object given a score function, list of loop segments to rebuild, and
	/// the fragment library to use.  The skip rate of each loop is _ignored_ (the
	/// decision of which loops to rebuild should be made outside the protocol).
	LoopRebuild(
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::loops::Loops Loops_in
	);

	~LoopRebuild();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief sets all the standard settings for LoopBuild
	void set_default_settings();

	/// @brief Get the protocol's Monte Carlo object
	protocols::moves::MonteCarloOP get_mc( core::pose::Pose & pose );

	/// @brief setting the mc object
	void set_default_mc( core::pose::Pose & pose );

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	/// @brief Are we extending the loops before running loop-rebuild?
	bool extended_loop();
	/// @brief Specify whether the protocol should extend the loops before running loop-rebuild
	void set_extended_loop( bool val );


	// ---------------- Acessors -------------------------------------------------

	/// @brief Were loops successfully built?
	bool get_success() { return success; }

	/// @brief
	bool get_random_loop_flag(){ return random_loop_flag_; }

	/// @brief
	void set_random_loop_flag(bool setting = true){ random_loop_flag_ = setting; }


	/// @brief
	int  get_allowed_failure_before_extend(){ return allowed_failure_before_extend_; }

	/// @brief
	void set_allowed_failure_before_extend(int setting ){ allowed_failure_before_extend_ = setting; }


	/// @brief
	int  get_allowed_failure_before_stop(){ return allowed_failure_before_stop_; }

	/// @brief
	void set_allowed_failure_before_stop(int setting ){ allowed_failure_before_stop_ = setting; }


	/// @brief
	bool get_abort_on_failed_loop(){ return abort_on_failed_loop_; }

	/// @brief
	void set_abort_on_failed_loop(bool setting = true){ abort_on_failed_loop_ = setting; }


	int allowed_failure_before_stop;
	/// @brief
	bool loop_model();

	/// @brief
	bool get_ccd_closure_exist();

	/// @brief
	bool get_desired_loops_exist();

	/// @brief
	core::Size desired_loops();

	/// @brief
	core::Real get_loop_combine_rate();

	/// @brief
	bool get_combine_if_fail_exist();

	/// @brief Gets the RMSd tolerence, the most a loop may move during reconstruction
	core::Real get_rmsd_tolerance();

	/// @brief Gets the chainbreak tolerence, the largest the chainbreak score may be to accept a loop conformation
	core::Real get_chain_break_tolerance();

private:


	/// @brief Rebuild all loops in the pose, chosen in a random order
	bool build_random_loops( core::pose::Pose & pose );

	/// @brief
	void set_looprlx_allow_move_map( int const & loop_begin, int const & loop_end, core::kinematics::MoveMap & mm );

	/// @brief Get the cutoff for score filtering loop conformations.  CURRENTLY UNUSED IN THE PROTOCOL
	core::Real get_score_filter_cutoff();

	/// @brief Select a loop at random from the input list
	bool select_one_loop(
		int nres,
		int & selected_loop,
		std::vector< int > & folded_loops,
		std::vector< int > & inter_res,
		int  & loop_begin,
		int  & loop_end,
		int  & cutpoint,
		bool & extend_this_loop,
		bool & are_loops_combined,
		int  & combine_interval,
		int  & loop_counter
	);


	/// @brief Extend the loop stems
	void barcode_extend_stems(
		core::pose::Pose & pose,
		int & barcst_extend_begin,
		int & barcst_extend_end,
		int & loop_begin,
		int & loop_end,
		int const & old_loop_begin,
		int const & old_loop_end,
		int const & nres,
		int const & selected_loop,
		int const & total_combine,
		int const & backward_combine
	);

	/// @brief If chainbreak score is greater than the tolerance, extend the loop
	void extend_barcode_regions_if_chain_break(
		core::pose::Pose & pose,
		int const & loop_begin,
		int const & loop_end,
		int & n_chain_break_fail,
		bool & is_chain_break,
		int & barcst_extend_begin, // output
		int & barcst_extend_end
	);

	/// @brief
	bool shorten_long_terminal_loop();


	/// @brief Is the RMS change between these two poses below the acceptible threshold?
	bool acceptable_rmsd_change( core::pose::Pose & pose1, core::pose::Pose & pose2 );

	/// @brief
	core::Real get_looprlx_cycle_ratio();

	/// @brief frag insertion + ccd close + minimize
	void build_loop_with_ccd_closure(
		core::pose::Pose & pose,
		int const & loop_begin,
		int const & loop_end,
		int & cutpoint,
		bool const & extend_this_loop
	);

	/// @brief Apply ccd moves to close the specified loop
	void fast_ccd_close_loops(
		core::pose::Pose & pose,
		int const & loop_begin,
		int const & loop_end,
		int const & cutpoint,
		core::kinematics::MoveMap & mm
	);

	// protocol-specific data
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::kinematics::MoveMapOP movemap_;
	protocols::moves::MonteCarloOP mc_;
	protocols::loops::Loops Loops_in_;
	std::vector< core::fragment::FragSetOP > frag_libs_;

	// parameters
	bool extend_loops;

	bool use_default_extend_loops;

	bool random_loop_flag_;

	bool only_remove_missing_density_flag_;

	int allowed_failure_before_extend_;

	int allowed_failure_before_stop_;

	bool abort_on_failed_loop_;

	core::Real m_Temperature_; // default temperature for monte carlo
	bool mc_created;

	bool success;
}; // class LoopRebuild

///////////////////////////////////////////////////////////////////
/// @brief class LoopRefine for fullatom loop refinement
///////////////////////////////////////////////////////////////////
class LoopRefine: public protocols::moves::Mover {
public:
	LoopRefine(
		protocols::loops::Loops Loops_in
	):Mover(),
		Loops_in_( Loops_in )
	{
		Mover::type("LoopRefine");
	}

	/// @brief Apply the loop-refine protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	protocols::loops::Loops Loops_in_;
};

} // protocols

#endif
