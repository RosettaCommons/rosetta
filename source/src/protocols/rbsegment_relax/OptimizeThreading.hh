// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_rbsegment_relax_optimize_threading_hh
#define INCLUDED_protocols_rbsegment_relax_optimize_threading_hh

#include <protocols/rbsegment_relax/OptimizeThreadingCreator.hh>

#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace rbsegment_relax {

class OptimizeThreadingMover : public moves::Mover {
public:
	OptimizeThreadingMover() : Mover(), nsteps_(1000), weight_(0.1) {
		scorefxn_ = new core::scoring::ScoreFunction();
		scorefxn_->set_weight( core::scoring::atom_pair_constraint , 1.0 );
	}

	virtual std::string get_name() const { return OptimizeThreadingMoverCreator::mover_name(); }
	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new OptimizeThreadingMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose );
	virtual void parse_my_tag(
			utility::tag::TagPtr const tag,
			moves::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Size nsteps_;
	core::Real weight_;
};


}
}

#endif
