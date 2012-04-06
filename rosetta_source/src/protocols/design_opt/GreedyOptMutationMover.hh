// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/design_opt/GreedyOptMutationMover.hh
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_design_opt_GreedyOptMutationMover_hh
#define INCLUDED_protocols_design_opt_GreedyOptMutationMover_hh
#include <protocols/design_opt/GreedyOptMutationMover.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#ifdef PYROSETTA
	#include <protocols/filters/Filter.hh>
#endif


namespace protocols {
namespace design_opt{

class GreedyOptMutationMover : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	GreedyOptMutationMover();
	//TODO: non-default ctor

	void apply( Pose & pose );
	protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new GreedyOptMutationMover ); }

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~GreedyOptMutationMover();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	protocols::filters::FilterOP filter() const;
	void filter( protocols::filters::FilterOP filter );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	bool dump_pdb() const;
	void dump_pdb( bool const dump_pdb );
	bool report_all() const;
	void report_all( bool const report_all );
	std::string sample_type() const;
	void sample_type( std::string const sample_type );
	void stopping_condition( protocols::filters::FilterOP f ){ stopping_condition_ = f; }
	protocols::filters::FilterOP stopping_condition() const{ return stopping_condition_; }
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::filters::FilterOP filter_;
	protocols::moves::MoverOP relax_mover_;
	std::string sample_type_;
	core::Real flip_sign_;
	bool report_all_;
	bool dump_pdb_;
	protocols::filters::FilterOP stopping_condition_; // dflt NULL ; if defined, stops greedy optimization when the filter's apply evaluates to true;
};


} // moves
} // protocols


#endif /*INCLUDED_protocols_design_opt_GreedyOptMutationMover_HH*/
