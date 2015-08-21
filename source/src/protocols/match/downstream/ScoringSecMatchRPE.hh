// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/GeometrySecMatchRPE.hh
/// @brief
/// @author Kui K. Chan, kuichan@u.washington.edu, Oct 09

#ifndef INCLUDED_protocols_match_downstream_ScoringSecMatchRPE_hh
#define INCLUDED_protocols_match_downstream_ScoringSecMatchRPE_hh


// Unit headers
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>
#include <protocols/match/downstream/ScoringSecMatchRPE.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {


class ScoringSecMatchRPE : public SecMatchResiduePairEvaluator
{

public:

	core::Size scoringEval_counter;

	//ScoringSecMatchRPE( std::string const & s_in, std::string const & pdb_file);

	//ScoringSecMatchRPE( std::string const & s_in);
	ScoringSecMatchRPE( std::string const & s_in, core::pose::Pose const & ref_pose);

	virtual ~ScoringSecMatchRPE();

	virtual
	core::Real
	getCutoff() const;

	core::scoring::ScoreFunctionCOP
	get_score_function() const;

	core::scoring::ScoreFunctionOP
	get_score_function();

	virtual
	void
	setPose( core::pose::Pose const & ref_pose );

	virtual
	bool
	evaluate_residues(
		core::conformation::Residue const & match_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	bool
	eval_cd_2b_residues(
		core::conformation::Residue const & match_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	bool
	eval_longRange_2b_residue(
		core::conformation::Residue const & match_res,
		core::conformation::Residue const & target_res
	) const;

	virtual
	bool
	require_all_target_residue_atom_coordinates() const;

	virtual
	bool
	require_target_atom_coordinate( Size target_atom_id ) const;

	virtual
	bool
	check2bsc(
		core::scoring::ScoreType sType,
		core::Real wts
	);


private:

	utility::vector1< core::scoring::ScoreType > secmatch_scotypes_cutoff_;
	utility::vector1< core::Real > secmatch_value_cutoff_;

	core::Real total_score_cutoff_;
	bool cutoff_flag_;
	bool cutoff_scoreType_flag_;

	core::pose::PoseOP ref_pose_;
	core::scoring::ScoreFunctionOP sfxn_;

	bool shortRange_;
	bool longRange_;
	core::pose::PoseOP cd_2b_pose_;


};


}
}
}

#endif
