// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.hh
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

#ifndef INCLUDED_protocols_forge_remodel_RemodelLoopMover_hh
#define INCLUDED_protocols_forge_remodel_RemodelLoopMover_hh

// unit headers
#include <protocols/forge/remodel/RemodelData.hh>
#include <protocols/forge/remodel/RemodelLoopMover.fwd.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>

// project headers
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>

// utility headers
#include <utility/vector1.hh>

#include <string>
#include <set>

#include <protocols/moves/MonteCarlo.fwd.hh>


namespace protocols {
namespace forge {
namespace remodel {


/// @brief Loop modeling protocol based on routines from Remodel and EpiGraft
///  packages in Rosetta++.
class RemodelLoopMover : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	typedef core::Real Real;
	typedef core::Size Size;
	typedef std::string String;

	typedef core::kinematics::MoveMap MoveMap;
	typedef core::fragment::FragSetCOP FragSetCOP;
	typedef core::fragment::FragSetOP FragSetOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

	typedef protocols::simple_moves::ClassicFragmentMoverOP ClassicFragmentMoverOP;
	typedef protocols::simple_moves::FragmentMoverOP FragmentMoverOP;
	typedef protocols::loops::Loop Loop;
	typedef protocols::loops::Loops Loops;
	typedef protocols::moves::MonteCarlo MonteCarlo;
	typedef protocols::moves::MoverOP MoverOP;

	typedef utility::vector1< FragSetOP > FragSetOPs;
	typedef utility::vector1< FragmentMoverOP > FragmentMoverOPs;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // construct/destruct


	/// @brief default constructor
	RemodelLoopMover();


	/// @brief loops constructor
	RemodelLoopMover( loops::LoopsOP const loops );


	/// @brief copy constructor
	RemodelLoopMover( RemodelLoopMover const & rval );


	/// @brief default destructor
	virtual
	~RemodelLoopMover();


public: // options


	void set_param_from_options();

	static void register_options();

private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	RemodelLoopMover & operator =( RemodelLoopMover const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;


	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;


public: // accessors


	/// @brief the ScoreFunction to use during modeling;
	ScoreFunction const & scorefunction() const;


	/// @brief get the false movemap
	/// @remarks All movemaps are generated with respect to this movemap.
	///  Any explicit False settings within this movemap will be retained in all
	///  descendant movemaps.
	inline
	MoveMap const & false_movemap() const {
		return false_movemap_;
	}


	/// @brief if linear chainbreak is <= this value, loop is considered closed
	///  (default 0.07)
	inline
	Real max_linear_chainbreak() const {
		return max_linear_chainbreak_;
	}


	/// @brief randomize loops prior to running main protocol? (default true)
	inline
	bool randomize_loops() const {
		return randomize_loops_;
	}


	/// @brief the allowed number of overall closure attempts before apply() exits
	///  (default 3)
	inline
	Size allowed_closure_attempts() const {
		return allowed_closure_attempts_;
	}

	/// @brief the number of loophash closure cycles to perform (default 8)
	inline
	Size loophash_cycles() const {
		return loophash_cycles_;
	}

	/// @brief the number of simultaneous closure cycles to perform (default 2)
	inline
	Size simultaneous_cycles() const {
		return simultaneous_cycles_;
	}


	/// @brief the number of independent closure cycles to perform (default 8)
	inline
	Size independent_cycles() const {
		return independent_cycles_;
	}


	/// @brief the maximum number of possible lockdown closure cycles to perform
	///  (default 30)
	inline
	Size boost_closure_cycles() const {
		return boost_closure_cycles_;
	}


	/// @brief the total number of "standard" (equal to simul + independent)
	///  to perform
	inline
	Size total_standard_cycles() const {
		return simultaneous_cycles_ + independent_cycles_ + loophash_cycles_;
	}

	/// @brief temperature for mc ( default 2.0 )
	inline
	Real temperature() const {
		return temperature_;
	}


public: // mutators

	// @brief for building repeat structures stored in private variable
	void repeat_generation_with_additional_residue(Pose & pose, Pose &repeat_pose);

	void
	repeat_sync(
		core::pose::Pose & repeat_pose,
		core::Size repeat_number
	);


	// @brief for updating repeat angles from a monomeric copy
	void repeat_propagation( Pose & pose, Pose & repeat_pose, Size repeat_number);


	/// @brief the ScoreFunction to use during modeling
	void scorefunction( ScoreFunction const & sfx );

	/// @brief adds the remodel data to thei remodelLoopMover object
	void remodelData( protocols::forge::remodel::RemodelData const remodel_data);

	/// @brief set the false movemap to use
	/// @remarks All movemaps are generated with respect to this movemap.
	///  Any explicit False settings within this movemap will be retained in all
	///  descendant movemaps.
	inline
	void false_movemap( MoveMap const & movemap ) {
		false_movemap_ = movemap;
	}


	/// @brief if linear chainbreak is <= this value, loop is considered closed
	inline
	void max_linear_chainbreak( Real const val ) {
		max_linear_chainbreak_ = val;
	}


	/// @brief randomize loops prior to running main protocol?
	inline
	void randomize_loops( bool const flag ) {
		randomize_loops_ = flag;
	}


	/// @brief the allowed number of overall closure attempts before apply() exits
	inline
	void allowed_closure_attempts( Size const attempts ) {
		allowed_closure_attempts_ = attempts;
	}

	/// @brief the number of loophash closure cycles to perform
	inline
	void loophash_cycles( Size const cycles ) {
		loophash_cycles_ = cycles;
	}

	/// @brief the number of simultaneous closure cycles to perform
	inline
	void simultaneous_cycles( Size const cycles ) {
		simultaneous_cycles_ = cycles;
	}


	/// @brief the number of independent closure cycles to perform
	inline
	void independent_cycles( Size const cycles ) {
		independent_cycles_ = cycles;
	}

	void
	set_user_provided_movers( utility::vector1< moves::MoverOP > const & movers );

	inline
	void user_provided_movers_apply_cycle( Size const cycle ) {
		user_provided_movers_apply_cycle_ = cycle; }


	/// @brief the maximum number of possible lockdown closure cycles to perform
	inline
	void boost_closure_cycles( Size const cycles ) {
		boost_closure_cycles_ = cycles;
	}

	/// @brief temperature for mc
	inline
	void temperature( Real const temp ) {
		temperature_ = temp;
	}

	inline
	void
	set_repeat_tail_length( Size const length ) {
		repeat_tail_length_ = length;
	}

	inline
	void
	set_keep_input_foldtree( bool const setting ){
		keep_input_foldtree_ = setting;
	}

	inline
	void
	set_seal_foldtree( bool const setting ){
		seal_foldtree_ = setting;
	}

public: // loop management


	/// @brief the loops to model
	inline
	loops::LoopsOP const loops() const {
		return loops_;
	}


	/// @brief set the loops to model
	inline
	void loops( loops::LoopsOP const loops ) {
		loops_ = loops;
	}


	/// @brief add a loop to model
	inline
	void add_loop( Loop const loop ) {
		loops_->add_loop( loop );
	}


public: // fragment management


	/// @brief add a fragment set
	void add_fragments( FragSetCOP fragset );


	/// @brief clear all fragment sets
	void clear_fragments();


public: // virtual main methods


	/// @brief apply defined moves to given Pose
	/// @remarks Sets protocols::moves::MS_SUCCESS upon successful closure of
	///  all loops, otherwise sets protocols::moves::FAIL_RETRY.
	virtual
	void apply( Pose & pose );

	virtual std::string get_name() const;

	//make this function public for experimental purpose

	/// @brief randomize loops
	void randomize_stage( Pose & pose );

protected: // loop modeling stages

	/// @brief find the smallest fragment size and insert a single such
	///  smallmer into each loop; for breaking up trapped trajectories
	/// @param[in,out] pose The pose to modify.
	/// @param[in] only_broken_loop If true, only insert into broken loops,
	///  otherwise insert into all. (default true)
	void insert_random_smallestmer_per_loop(
		Pose & pose,
		bool const only_broken_loops = true
	);

	void set_segment_stage(
		Pose & pose,
		MonteCarlo & mc
	);


	void loophash_stage(
		Pose & pose,
		MonteCarlo & mc,
		Real const cbreak_increment
	);
	/// @brief abinitioo stage: better control of fragments. No loop closure tried
	void abinitio_stage(
		Pose & pose,
		Size const fragmentSize,
		MoveMap const movemap,
		ScoreFunctionOP sfxOP,
		Size const max_outer_cycles,
		Size const max_inner_cycles,
		std::set<Size> const & disallowedPos,
		bool const recover_low,
		std::string stage_name,
		bool const smoothMoves,
		Real const fragScoreThreshold
	);
	/// @brief quick trip to FA to push parts of the molecular apart to a reasonable distance. I'm noticing an overcollapse in centroid
	void fa_relax_stage(
		Pose & pose
	);

	/// @ brief sets up ncs constraints
	protocols::simple_moves::symmetry::SetupNCSMover generate_ncs_csts(Pose & pose);

	/// @brief small_move stage: only small moves
	void small_move_stage(
		Pose & pose,
		MoveMap const movemap,
		ScoreFunctionOP sfxOP,
		Size const max_outer_cycles,
		Size const max_inner_cycles,
		bool const recover_low,
		Real const h_range,
		Real const e_range,
		Real const l_range);


	/// @brief simultaneous stage: multiple loop movement prior to MC accept/reject
	void simultaneous_stage(
		Pose & pose,
		MonteCarlo & mc,
		Real const cbreak_increment
	);


	/// @brief independent stage: single loop movement prior to MC accept/reject
	void independent_stage(
		Pose & pose,
		MonteCarlo & mc,
		Real const cbreak_increment
	);


	/// @brief lockdown stage: close loops within some threshold
	///  w/ smallest-mer (typically 1-mer) + ccd_move only
	void boost_closure_stage(
		Pose & pose,
		MonteCarlo & mc,
		Real const cbreak_increment
	);


protected: // loops


	/// @brief determine which loops need modeling wrt to given Pose
	/// @remarks Skips closed loops and shuffles the order of the remaining
	///  loops.
	loops::LoopsOP determine_loops_to_model( Pose & pose );


	/// @brief check all loops for closure criteria
	/// @param[in] pose The pose being checked.
	/// @param[in] show_in_tracer Output state of each loop to tracer?
	/// @return true if all criteria pass, false otherwise
	bool check_closure_criteria(
		Pose & pose,
		bool const show_in_tracer = false
	);


protected: // fragments


	/// @brief return fragment movers for the list of internally kept fragment sets,
	///  1 fragment mover for each fragment set
	/// @param[in] movemap Use this movemap when initializing fragment movers.
	/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
	///  size is this number.  If zero, uses all fragment sets.
	FragmentMoverOPs create_fragment_movers(
		MoveMap const & movemap,
		Size const largest_frag_size = 0
	);
	/// @brief same as create_fragment_movers except with set size
	FragmentMoverOPs create_fragment_movers_limit_size(
		MoveMap const & movemap,
		Size const frag_size,
		std::set<Size> const & disallowedPos,
		bool const smoothMoves,
		Real const fragScoreThreshold = 999.0
	);

	/// @brief append fragment movers for the list of internally kept fragment sets,
	///  1 fragment mover for each fragment set
	/// @param[in] movemap Use this movemap/s when initializing fragment movers.
	/// @param[out] frag_movers Append fragment movers to this list.
	/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
	///  size is this number.  If zero, uses all fragment sets.
	void create_fragment_movers(
		MoveMap const & movemap,
		FragmentMoverOPs & frag_movers,
		Size const largest_frag_size = 0
	);

	/// @brief create per-loop fragment movers: 1 fragment mover for each loop (uses
	///  movemaps to lock down non-loop residues)
	/// @param[in] loops The loops to use.
	/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
	///  size is this number.  If zero, uses all fragment sets.
	FragmentMoverOPs create_per_loop_fragment_movers(
		loops::LoopsOP const loops,
		Size const largest_frag_size = 0
	);


protected: // movemap


	/// @brief enforce settings in the false movemap
	void enforce_false_movemap( MoveMap & movemap );


	/// @brief mark bb/chi torsions of multiple loops moveable in a movemap
	/// @param[in] loops The loops to use.
	/// @param[out] movemap The movemap to modify.
	/// @param[in] allow_omega Allow bb omega to move? (should be yes when
	///  doing either fragment insertion or scoring function has omega
	///  tether, otherwise should probably be no)
	void mark_loops_moveable(
		loops::LoopsOP const loops,
		MoveMap & movemap,
		bool const allow_omega
	);


	/// @brief mark bb/chi torsion of a single loop moveable in movemap
	/// @param[in] loops The loop to use.
	/// @param[out] movemap The movemap to modify.
	/// @param[in] allow_omega Allow bb omega to move? (should be yes when
	///  doing either fragment insertion or scoring function has omega
	///  tether, otherwise should probably be no)
	void mark_loop_moveable(
		Loop const & loop,
		MoveMap & movemap,
		bool const allow_omega
	);


	/// @brief count number of residues with moveable backbone torsions in the
	///  given range [left, right]
	Size count_moveable_residues(
		MoveMap const & movemap,
		Size const left,
		Size const right
	);
private:
	/// @brief initializes pose with starting pose
	void set_starting_pdb(Pose & pose);

	/// @brief initializes pose with starting sequence
	void set_starting_sequence(Pose & pose);

	/// @brief initializes pose with ideal helices.
	void set_ideal_helices(Pose & pose);

	/// @setup allowed positions per stage
	std::set<core::Size> generate_residues_to_sample(bool chooseSubsetResidues, Pose & pose, Size fragmentSize);

private: // data

	/// @brief the score function to use
	ScoreFunctionOP sfx_;

	/// @remodelData, used when constraints are defined through Remodel
	//So you can get at everything in the blueprint data
	protocols::forge::remodel::RemodelData remodel_data_;

	/// @brief the false movemap to use
	/// @remarks All movemaps are generated with respect to this movemap.
	///  Any False explicit settings within this movemap will be retained in all
	///  descendant movemaps.
	MoveMap false_movemap_;


	/// @brief list of loops to model
	loops::LoopsOP loops_;


	/// @brief if linear chainbreak is <= this value, loop is considered closed
	///  (default 0.07)
	Real max_linear_chainbreak_;


	/// @brief randomize loops prior to running main protocol? (default true)
	bool randomize_loops_;

	Size repeat_tail_length_;


	/// @brief the allowed number of overall closure attempts before apply() exits
	///  (default 3)
	Size allowed_closure_attempts_;

	/// @brief the number of loophash closure cycles to perform (default 8)
	Size loophash_cycles_;

	/// @brief the number of simultaneous closure cycles to perform (default 2)
	Size simultaneous_cycles_;


	/// @brief the number of independent closure cycles to perform (default 8)
	Size independent_cycles_;

	/// @brief collection of movers that mess with the pose
	/// every n fragment/ccd steps during inneriterations of simultaneous,
	/// and independent stage, where n is determined by the variable below
	/// gives the user a chance to add their own flavor to the fragment sampling
	/// done by this protocol
	utility::vector1< moves::MoverOP > user_provided_movers_;

	/// @brief determines how often the above movers get called
	/// during fragment insertion/ccd steps (default 3 )
	Size user_provided_movers_apply_cycle_;

	/// @brief the maximum number of possible boost closure cycles to perform
	///  (default 30)
	Size boost_closure_cycles_;

	/// @brief temperature for mc
	Real temperature_;


	/// @brief fragment sets to use
	FragSetOPs fragsets_;

	/// @brief local copy of repeat pose
	Pose repeat_pose_;

	/// @brief use this for repeat construction in conjunction with symmetry
	Size unit_length_;
	Size repeat_length_;


	/// @brief RemodelGlobalFrame for helical definitions
	RemodelGlobalFrame RGF_;

	/// @brief whether to keep the FoldTree untouched, i.e. trust the
	///user to supply a decent FoldTree
	bool keep_input_foldtree_;

	/// @brief whether to (true) generate a sealed fold tree or (false) keep the cut fold tree (default=true)
	bool seal_foldtree_;
};


} // remodel
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_remodel_RemodelLoopMover_HH */
