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
#include <protocols/protein_interface_design/filters/DeltaFilter.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#include <protocols/filters/Filter.hh>


namespace protocols {
namespace design_opt{

class GreedyOptMutationMover : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	GreedyOptMutationMover();
	GreedyOptMutationMover(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		protocols::filters::FilterOP filter,
		core::Real filter_delta = 0,
		std::string sample_type = "low",
		bool stop_before_condition = false,
		bool skip_best_check = false,
		bool dump_pdb = false,
		bool dump_table = false,
		bool rtmin = false,
		bool parallel = false,
		core::Size diversify_lvl = core::Size( 1 ),
		protocols::filters::FilterOP stopping_condition = protocols::filters::FilterOP( NULL )
	);

	bool pose_coords_are_same( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
	void clear_cached_data();
	void apply( Pose & pose );
	protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new GreedyOptMutationMover ); }

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &,
			protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~GreedyOptMutationMover();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	core::Real filter_delta() const;
	void filter_delta( core::Real filter_delta );
	protocols::filters::FilterOP filter() const;
	void filter( protocols::filters::FilterOP filter );
	bool stop_before_condition() const;
	void stop_before_condition( bool const stop_before_condition );
	bool skip_best_check() const;
	void skip_best_check( bool const skip_best_check );
	bool dump_pdb() const;
	void dump_pdb( bool const dump_pdb );
	bool dump_table() const;
	void dump_table( bool const dump_table );
	void dump_scoring_table( std::string filename, core::pose::Pose const & ref_pose  ) const;
	std::string sample_type() const;
	void sample_type( std::string const sample_type );
	core::Size diversify_lvl() const;
	void diversify_lvl( core::Size const diversify_lvl );
	void stopping_condition( protocols::filters::FilterOP f ){ stopping_condition_ = f; }
	protocols::filters::FilterOP stopping_condition() const{ return stopping_condition_; }
	bool rtmin() const;
	void rtmin( bool const b );
	bool parallel() const;
	void parallel( bool const b );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::moves::MoverOP relax_mover_;
	protocols::filters::FilterOP filter_;
	core::Real filter_delta_;
	utility::vector1<protocols::protein_interface_design::filters::DeltaFilterOP> reset_delta_filters_;
	std::string sample_type_;
	bool stop_before_condition_;
	bool skip_best_check_;
	bool dump_pdb_;
	bool dump_table_;
	core::Size diversify_lvl_;
	core::Real flip_sign_;
	protocols::filters::FilterOP stopping_condition_; // dflt NULL ; if defined, stops greedy optimization when the filter's apply evaluates to true;
	utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, core::Real > > > > seqpos_aa_val_vec_;
	core::pose::Pose ref_pose_;
	bool rtmin_; //dflt false; should we rtmin after packing?
	bool parallel_; //parallelize pointmut calc with MPI?
	core::Real design_shell_;//dflt -1 to only allow pointmutations, higher allows suroundings to be designed as well
	core::Real repack_shell_;
};


} // design_opt
} // protocols


#endif /*INCLUDED_protocols_design_opt_GreedyOptMutationMover_HH*/
