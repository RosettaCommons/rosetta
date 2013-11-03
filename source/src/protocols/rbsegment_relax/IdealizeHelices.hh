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


#ifndef INCLUDED_protocols_rbsegment_relax_idealize_helices_hh
#define INCLUDED_protocols_rbsegment_relax_idealize_helices_hh

#include <protocols/rbsegment_relax/IdealizeHelicesCreator.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>


namespace protocols {
namespace rbsegment_relax {

class IdealizeHelicesMover : public moves::Mover {
public:
	IdealizeHelicesMover() : Mover(), cst_weight_(1.0), cst_width_(0.0) {
		scorefxn_ = new core::scoring::ScoreFunction();
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
		scorefxn_->set_weight( core::scoring::coordinate_constraint , cst_weight_ );
	}

	virtual std::string get_name() const { return IdealizeHelicesMoverCreator::mover_name(); }
	moves::MoverOP clone() const { return( protocols::moves::MoverOP( new IdealizeHelicesMover( *this ) ) ); }

	virtual void apply( core::pose::Pose & pose );

	virtual void parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap &data,
			filters::Filters_map const &filters,
			moves::Movers_map const &movers,
			core::pose::Pose const & pose );

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real cst_weight_,cst_width_;
	utility::vector1< std::pair<int,int> > helices_;
};


}
}

#endif
