// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_ScoreEvaluator_hh
#define INCLUDED_protocols_simple_filters_ScoreEvaluator_hh


// Unit Headers
// #include <protocols/simple_filters/ScoreEvaluator.fwd.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>
// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/rms_util.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_filters/RDC_Evaluator.fwd.hh>
#include <utility/vector1.hh>


// ObjexxFCL Headers

// Utility headers


//// C++ headers

namespace protocols {
namespace simple_filters {

/// @brief that rewrites the whole pss struct all previous entries will be lost... probably not what one wants...
class ScoreEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	ScoreEvaluator( std::string tag, core::scoring::ScoreFunctionOP scorefxn, bool fullname=false );
	~ScoreEvaluator();
	virtual core::Real
	apply( core::pose::Pose& pose ) const;
	virtual bool applicable( core::pose::Pose const&pose ) const;
protected:
	core::scoring::ScoreFunctionOP scorefxn_;

};

class TruncatedScoreEvaluator : public ScoreEvaluator {
public:

	TruncatedScoreEvaluator( std::string tag, core::scoring::ResidueSelectionVector const&, core::scoring::ScoreFunctionOP scorefxn=NULL, bool fullname=false );
	virtual ~TruncatedScoreEvaluator();
	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;
private:
	core::scoring::ResidueSelectionVector selection_;
	mutable core::scoring::ResidueSelectionVector exclude_list_;
	mutable Size nres_;
	SelectRDC_EvaluatorOP rdcs_;
};

}
}

#endif
