// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_rbsegment_relax_optimize_threading_hh
#define INCLUDED_protocols_rbsegment_relax_optimize_threading_hh


#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>


namespace protocols {
namespace rbsegment_relax {

class OptimizeThreadingMover : public moves::Mover {
public:
	OptimizeThreadingMover() : Mover(),
		nsteps_(1000), rebuild_cycles_(200), max_shift_(4), weight_(0.1), temperature_(2.0), recover_low_(true), step_penalty_(false), greedy_(false), native_(/* NULL */) {
		scorefxn_ = utility::pointer::make_shared< core::scoring::ScoreFunction >();
		scorefxn_->set_weight( core::scoring::atom_pair_constraint , 1.0 );
		scorefxn_sampling_ = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth");
		if ( scorefxn_sampling_->get_weight(core::scoring::atom_pair_constraint) == 0 ) {
			scorefxn_sampling_->set_weight( core::scoring::atom_pair_constraint , 1.0 );
		}
		if ( scorefxn_sampling_->get_weight(core::scoring::dihedral_constraint) == 0 ) {
			scorefxn_sampling_->set_weight( core::scoring::dihedral_constraint , 1.0 );
		}

		loops_ = utility::pointer::make_shared< protocols::loops::Loops >();
	}

	moves::MoverOP clone() const override { return( utility::pointer::make_shared< OptimizeThreadingMover >( *this ) ); }

	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	// helper function
	void rebuild_unaligned(core::pose::Pose &pose);

	core::scoring::ScoreFunctionOP scorefxn_, scorefxn_sampling_;
	core::Size nsteps_, rebuild_cycles_, max_shift_;
	core::Real weight_,temperature_;
	bool recover_low_,step_penalty_,greedy_;
	utility::vector1< std::string > chains_;

	protocols::loops::LoopsOP loops_;

	core::pose::PoseOP native_; // "cheating" ... optimize RMS to native
};


}
}

#endif
