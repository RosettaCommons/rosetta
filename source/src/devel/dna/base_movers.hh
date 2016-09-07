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

#ifndef INCLUDED_devel_dna_base_movers_hh
#define INCLUDED_devel_dna_base_movers_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <devel/cartesian_frags/DNA_FragLib.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility>
#include <utility/pointer/owning_ptr.hh>

#include <utility/vector1.hh>


namespace devel {
namespace dna {

///////////////////////////////////////////////////////////////////////////////
/// @brief  A mover that makes base-pair fragment insertions
/// @details  Uses a cartesian-fragment library for the fragments to insert and also for the
/// ensuing patch-up process where the chainbreaks are fixed up.


class BasePairMover : public protocols::moves::Mover {

public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef cartesian_frags::DNA_FragLibOP DNA_FragLibOP;

public:
	BasePairMover(
		DNA_FragLibOP lib,
		Real const frag_dev_threshold,
		Size const max_tries,
		Real const max_score_increase,
		ScoreFunctionOP scorefxn
	):
		lib_(std::move( lib )),
		frag_dev_threshold_( frag_dev_threshold ),
		max_tries_( max_tries ),
		max_score_increase_( max_score_increase ),
		scorefxn_(std::move( scorefxn ))
	{
		protocols::moves::Mover::type( "BasePairMover" );
	}


	void
	apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	DNA_FragLibOP lib_;
	Real frag_dev_threshold_;
	Size max_tries_;
	Real max_score_increase_;
	ScoreFunctionOP scorefxn_;

};

typedef utility::pointer::shared_ptr< BasePairMover > BasePairMoverOP;


///////////////////////////////////////////////////////////////////////////////
/// @brief  A mover that makes base-step fragment insertions
/// @details  Uses a cartesian-fragment library for the fragments to insert and also for the
/// ensuing patch-up process where the chainbreaks are fixed up.

class BaseStepMover : public protocols::moves::Mover {

public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef cartesian_frags::DNA_FragLibOP DNA_FragLibOP;

public:
	BaseStepMover(
		DNA_FragLibOP lib,
		Real const frag_dev_threshold,
		Size const max_tries,
		Real const max_score_increase,
		ScoreFunctionOP scorefxn
	):
		lib_(std::move( lib )),
		frag_dev_threshold_( frag_dev_threshold ),
		max_tries_( max_tries ),
		max_score_increase_( max_score_increase ),
		scorefxn_(std::move( scorefxn ))
	{
		protocols::moves::Mover::type( "BaseStepMover" );
	}


	void
	apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:
	DNA_FragLibOP lib_;
	Real frag_dev_threshold_;
	Size max_tries_;
	Real max_score_increase_;
	ScoreFunctionOP scorefxn_;

};

typedef utility::pointer::shared_ptr< BaseStepMover > BaseStepMoverOP;

} // ns dna
} // ns devel

#endif
