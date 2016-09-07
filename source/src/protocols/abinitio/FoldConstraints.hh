// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
///
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_FoldConstraints_hh
#define INCLUDED_protocols_abinitio_FoldConstraints_hh


// Unit Headers
#include <protocols/abinitio/FoldConstraints.fwd.hh>

// Package Headers
#include <protocols/simple_moves/FragmentMover.fwd.hh>
//#include <protocols/simple_moves/GunnCost.fwd.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
//#include <protocols/abinitio/ConstraintFragmentMover.fwd.hh>
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <core/kinematics/FoldTree.fwd.hh>

//#include <core/scoring/EnergyMap.fwd.hh>

#include <protocols/moves/Mover.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>

#include <protocols/constraints_additional/MaxSeqSepConstraintSet.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

/// @brief extension of ClassicAbinitio Protocol to adapt the folding process for the presence of distance constraints
/// @detail Main Function: switch distance constraints based on distance in the FoldTree ( i.e., in sequence for simple FoldTrees )
///     This is achieved by replacing the pose's ConstraintSet with the special purpose class MaxSeqSepConstraintSet
///     the latter class will only score constraints that are sufficiently close in FoldTree/Sequence
///     ( as controlled by the threshold with set_max_seq_sep()  )
///     the protocol ranks up the max_seq_sep parameter while folding proceeds through the stages.
///     to this extend it overloads methods prepare_stageX() do_stage1_cycles()
///
///     the other substantial difference to ClassicAbinitio is that minimizations are carried out.
///     method min_trial() is called each time the max_seq_sep is changed. ( inhibit: -no_minimize )
class FoldConstraints : public ClassicAbinitio {
public:
	typedef ClassicAbinitio Parent;

public:
	/// @brief c'stor from Movers
	FoldConstraints(
		simple_moves::FragmentMoverOP brute_move_small,
		simple_moves::FragmentMoverOP brute_move_large,
		simple_moves::FragmentMoverOP smooth_move_small,
		int dummy /* otherwise the two constructors are ambigous */
	);

	/// @brief c'stor from FragSets --- ClassicFragmentMover and SmoothFragmentMover will be created
	FoldConstraints(
		core::fragment::FragSetCOP fragset3mer,
		core::fragment::FragSetCOP fragset9mer,
		core::kinematics::MoveMapCOP movemap
	);

	/// @brief Explicit copy constructor to handle OPs.
	FoldConstraints( FoldConstraints const & src );

	/// @brief Explicit destructor to handle OPs
	~FoldConstraints() override;

	/// @brief ...
	moves::MoverOP clone() const override;

	/// @brief run the protocol
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

	/// @brief sets the usual scores ( score0,score1, score2/5 etc. ) and additionally atom_pair_constraints to 1.0
	void set_default_scores() override;

	/// @brief
	void set_constraint_weight( core::Real setting ) {
		constraint_weight_ = setting;
		set_score_weight( core::scoring::atom_pair_constraint, constraint_weight_ );
	};

	void set_default_options() override;
	static void register_options();

	void set_show_viol_level( core::Size setting ) {
		show_viol_level_ = setting;
	}
protected:
	//overload some methods of ClassicAbinitio to change the MaxSeqSep of the Constraints throughout the protocol
	bool prepare_stage1( core::pose::Pose& pose ) override;
	bool prepare_stage2( core::pose::Pose& pose ) override;
	bool prepare_stage4( core::pose::Pose& pose ) override;
	bool prepare_loop_in_stage3( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) override;
	bool prepare_loop_in_stage4( core::pose::Pose &pose, Size loop_iteration, Size total_iterations ) override;

	bool do_stage1_cycles( core::pose::Pose& pose ) override;
	bool do_stage2_cycles( core::pose::Pose& pose ) override;

	virtual void setup_default_min_move();

	//@brief change the movemap ( is propagated to mover-objects )
	void set_movemap ( core::kinematics::MoveMapCOP mm ) override;

	void set_min_move( protocols::simple_moves::MinMoverOP mm);

	protocols::simple_moves::MinMover& min_move() {
		return *min_move_;
	}

	void min_trial( core::pose::Pose& pose );

	virtual void set_max_seq_sep( core::pose::Pose& pose, Size setting );

	core::Real max_seq_sep_fudge() const {
		return max_seq_sep_fudge_;
	}

	void max_seq_sep_fudge( core::Real setting ) {
		max_seq_sep_fudge_ = setting;
	}

	constraints_additional::MaxSeqSepConstraintSet const& constraints() {
		return *constraints_;
	}

	void
	bIgnoreSequenceSeparation( bool setting ) {
		bIgnoreSequenceSeparation_ = setting;
	}

	bool
	bIgnoreSequenceSeparation() {
		return bIgnoreSequenceSeparation_;
	}

	Size total_res( core::pose::Pose const& pose ) const;

	void set_seq_sep_stage1 ( core::Real setting ) {
		seq_sep_stage1_ = setting;
	}

	void set_seq_sep_stage3 ( core::Real setting ) {
		seq_sep_stage3_ = setting;
	}

private:
	core::Real
	evaluate_constraint_energy( core::pose::Pose& pose, core::scoring::ScoreFunction const& ) const;

	protocols::simple_moves::MinMoverOP min_move_;

	constraints_additional::MaxSeqSepConstraintSetOP constraints_;
	core::Real constraint_weight_;

	bool bMinTrial_;
	bool bIgnoreSequenceSeparation_;
	Size run_;

	core::Real max_seq_sep_fudge_;
	core::Real seq_sep_stage1_;
	core::Real seq_sep_stage3_;
	core::Real seq_sep_stage4_;

	core::Real start_ramp_cstweight_;
	core::Size ramp_cst_cycles_;
	core::Size ramp_iterations_;
	bool bSkipOnNoViolation_;

	//just for screen output: how verbose should it be
	Size show_viol_level_;

	//abolish run in stage2 if constraint threshold is violated -- '0' = inactive
	Size constraint_threshold_;
};


} //abinitio
} // protocols

#endif
