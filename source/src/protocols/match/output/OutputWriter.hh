// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/OutputWriter.fwd.hh
/// @brief  Forward declaration of class to write output matches that have (presumably) passed the output filters.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_OutputWriter_hh
#define INCLUDED_protocols_match_output_OutputWriter_hh

// Unit headers
#include <protocols/match/output/OutputWriter.fwd.hh>

// Package Headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/MatcherTask.fwd.hh>
#include <protocols/match/output/MatchEvaluator.fwd.hh>
#include <protocols/match/output/MatchScoreWriter.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ unit headers
#include <map>

namespace protocols {
namespace match {
namespace output {

class OutputWriter : public utility::pointer::ReferenceCount
{
public:
	OutputWriter();

	virtual ~OutputWriter();

	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP
	cst_io() const;

	virtual
	void
	record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer ) = 0;

	/// @brief evaluator and score writer are not passed in because single-downstream-position match
	///currently have no way of being evaluated
	virtual
	void
	record_match( match_dspos1 const & m ) = 0;

	/// @brief determine if any upstream res are at the same scaffold
	/// position, i.e. if one of them is a backbone interaction
	/// the redundant_upstream_res map is a mapping from the redundant
	/// geometric constraint id of the redundant residue to the geomcst_id of
	/// the "nonredundant" res i.e. if cstres 1 happens to be a cys at
	/// position 10 and cstres 3 is a gly at position 10 (and it's backbone is
	/// used in satisfying constraint 3), the mapping will contain the 3,1 pair.
	void
	determine_redundant_upstream_matchres(
		match_dspos1 const & m,
		std::map< core::Size, core::Size > & redundant_upstream_res
	) const;

	virtual
	void
	prepare_for_output_writing(){}

	virtual
	void
	end_output_writing(){}

	virtual
	void
	initialize_from_matcher_task(
		MatcherTaskCOP mtask
	);


private:

	protocols::toolbox::match_enzdes_util::EnzConstraintIOCOP cst_io_;

};

}
}
}

#endif
