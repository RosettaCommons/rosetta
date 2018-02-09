// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LocalInserter.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_loophash_LocalInserter_hh
#define INCLUDED_protocols_loophash_LocalInserter_hh

//Unit
#include <protocols/loophash/LocalInserter.fwd.hh>

//Core
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/MinimizerOptions.hh>

//Protocols
#include <protocols/loophash/BackboneDB.hh>

//Utility
#include <utility/pointer/ReferenceCount.hh>

//C++
#include <string>
#include <vector>

#include <utility/vector1.hh>

namespace protocols {
namespace loophash {

/// @brief Manages the insertion of an arbitrary length of backbone in a local 
/// manner.
/// @details This is a pure virtual superclass, and the intention is that 
/// different subclasses can use different methods to keep the insertion local.
class LocalInserter : public utility::pointer::ReferenceCount {
public:

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~LocalInserter() override;

	LocalInserter(){
	}

	/// @brief Insert a backbone segment.
	/// @param[out] start_pose
	/// @param[in] original_pose
	/// @param[in] new_bs The backbone segment to insert
	/// @param[in] respos The residue where the start of the given backbone 
	/// segment should be inserted.
	virtual
	core::Real
	make_local_bb_change(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) = 0;

	/// @brief Insert a backbone segment, and close any chainbreaks in the region 
	/// where the segment is being inserted.
	/// @param[out] start_pose
	/// @param[in] original_pose
	/// @param[in] new_bs The backbone segment to insert
	/// @param[in] respos The residue where the start of the given backbone 
	/// segment should be inserted.
	virtual
	core::Real
	make_local_bb_change_close_gaps(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) = 0;

	/// @brief Closes many gaps outside of ir and jr.
	/// @details Will die if gap exists between ir and jr.
	/// @param[out] start_pose
	/// @param[in] original_pose
	/// @param[in] new_bs The backbone segment to insert
	/// @param[in] respos The residue where the start of the given backbone 
	/// segment should be inserted.
	virtual
	core::Real
	make_local_bb_change_include_cut(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) = 0;

private:

};

/// @brief Insert a backbone segment and use minimization with coordinate 
/// constraints to keep the insertion local.
class LocalInserter_SimpleMin : public LocalInserter{
public:
	LocalInserter_SimpleMin():
		LocalInserter(),
		scorefxn_rama_cst_( core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
		options_( "dfpmin", 0.2, true , false ),
		scorefxn_cen_cst_( core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
		options2_( "dfpmin", 0.02,true , false )
	{
		set_default_score_functions();
	}

	/// @brief Set the score function for the first round of minimization
	/// during a loophash insert.
	void
	scorefxn_rama_cst(
		core::scoring::ScoreFunction const &  scorefxn
	);

	/// @brief Set the score function for the second round of minimization
	/// during a loophash insert.
	void
	scorefxn_cen_cst(
		core::scoring::ScoreFunction const & scorefxn
	);

	/// @details This method uses the following algorithm:
	/// - Constrain atoms outside the region where the insertion will occur to 
	///   their initial coordinates.
	/// - Apply the torsions from the given backbone segment.
	/// - Minimize first with only the rama term
	/// - Minimize again with a hard-coded centroid-like score function (see 
	///   LocalInserter_SimpleMin::set_default_score_functions() for details) and 
	///   a stricter convergence tolerance.
	core::Real
	make_local_bb_change(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) override;

	/// @details Similar to make_local_bb_change(), but before the torsions from 
	/// the given backbone segment are inserted, the bond lengths and angles for 
	/// those residues are idealized, so that if there was a chainbreak, it would 
	/// be fixed.
	/// @note This whole function, with the exception of the "idx" for-loop, is 
	/// duplicated from make_local_bb_change().  This could be easily refactored.
	core::Real
	make_local_bb_change_close_gaps(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) override;

	/// @details It looks to me like this method applies coordinate constraints 
	/// and inserts the given backbone segment, but never does any minimization.  
	core::Real
	make_local_bb_change_include_cut(
		core::pose::Pose &start_pose,
		const core::pose::Pose &original_pose,
		const protocols::loophash::BackboneSegment &new_bs,
		core::Size res_pos
	) override;

private:
	/// @brief Define hard-coded default scorefunctions for the minimization 
	/// steps.
	void set_default_score_functions();

	// the scorefunctions themselves
	core::scoring::ScoreFunctionOP scorefxn_rama_cst_;
	core::optimization::MinimizerOptions options_;

	core::scoring::ScoreFunctionOP scorefxn_cen_cst_;
	core::optimization::MinimizerOptions options2_;
};


}

}

#endif

