// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/OtherContextScoreFunction.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Florian Richter (floric@u.washington.edu)

#ifndef INCLUDED_protocols_flexpack_OtherContextScoreFunction_hh
#define INCLUDED_protocols_flexpack_OtherContextScoreFunction_hh

/// Unit headers
/// 1. CREATE THE FORWARD DECLARATION FILE
/// 2. #INCLUDE IT HERE

/// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>  // GET OUT OF THE HABIT OF #INCLUDING POSE.HH IN HEADER FILES

#include <utility/vector1.hh>



namespace protocols {
namespace flexpack {

class OtherContextScoreFunction : public core::scoring::ScoreFunction
{

public:
	typedef core::scoring::ScoreFunction parent;

public:
	OtherContextScoreFunction();
	~OtherContextScoreFunction();/// ALWAYS EXPLICITLY DEFINE A DSTOR -- *ESPECIALLY* IF THE CLASS CONTAINS AN OP TO SOME OTHER CLASS

	OtherContextScoreFunction(
		core::pose::Pose const & context_pose
	);

private:
	// private assignment and copy constructors to avoid discarding subtype information

	OtherContextScoreFunction( OtherContextScoreFunction const & );

	OtherContextScoreFunction &
	operator=( OtherContextScoreFunction const & );

public:

	// TODO: This class probably needs to have clone() and assign() methods if it's going to see non-trivial use.

	void
	set_context_pose( core::pose::Pose const & pose );

	void pre_scoring(); //call this before scoring.

	void
	eval_cd_1b(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


	/// @brief accumulate unweighted interaction energies between rsd1 and rsd2
	/// for all short ranged context dependent two body energies contained in this scorefunction
	void
	eval_cd_2b(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2,
		core::pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


	void
	eval_cd_intrares_energy(
		core::conformation::Residue const & rsd,
		core::pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const;


private:

	core::pose::PoseOP context_pose_;
	bool scored_context_pose_;
};

/// CREATE A .FWD.HH FILE AND PUT THIS TYPEDEF THERE
typedef utility::pointer::owning_ptr< OtherContextScoreFunction > OtherContextScoreFunctionOP;
typedef utility::pointer::owning_ptr< OtherContextScoreFunction const > OtherContextScoreFunctionCOP;

}
}

#endif
