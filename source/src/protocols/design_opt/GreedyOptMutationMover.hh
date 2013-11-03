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
#include <protocols/simple_filters/DeltaFilter.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
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
		utility::vector1< protocols::filters::FilterOP > filters,
		utility::vector1< std::string > sample_types,
		utility::vector1< core::Real > filter_deltas,
		bool dump_pdb = false,
		bool dump_table = false,
		bool parallel = false,
    bool stop_before_condition = false,
    bool skip_best_check = false,
    bool rtmin = false,
    bool shuffle_order = false,
    bool diversify = false,
    bool incl_nonopt = false,
		protocols::filters::FilterOP stopping_condition = protocols::filters::FilterOP( NULL )
	);

	bool pose_coords_are_same( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
	void filter_seqpos_pareto_opt_ptmuts();
	void filter_pareto_opt_poses();
	void clear_cached_data();
	void calc_pfront_poses_filter_ranks();
	void dump_scoring_table( std::string filename, core::pose::Pose const & ref_pose ) const;
	void apply( Pose & pose );
	protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new GreedyOptMutationMover ); }

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void add_filter(
		protocols::filters::FilterOP,
		std::string const sample_type,
		core::Real filter_delta
	);
	virtual ~GreedyOptMutationMover();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	utility::vector1< protocols::filters::FilterOP > filters() const;
	void filters( utility::vector1< protocols::filters::FilterOP > filters );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	bool dump_pdb() const;
	void dump_pdb( bool const dump_pdb );
	bool dump_table() const;
	void dump_table( bool const dump_table );
	bool parallel() const;
	void parallel( bool const parallel );
	utility::vector1< std::string > sample_types() const;
	void sample_types( utility::vector1< std::string > const sample_types );
	utility::vector1< core::Real > filter_deltas() const;
	void filter_deltas( utility::vector1< core::Real > const filter_deltas );
	void stopping_condition( protocols::filters::FilterOP f ){ stopping_condition_ = f; }
	protocols::filters::FilterOP stopping_condition() const{ return stopping_condition_; }
  bool stop_before_condition() const;
  void stop_before_condition( bool const stop_before_condition );
  bool skip_best_check() const;
  void skip_best_check( bool const skip_best_check );
  utility::vector1< protocols::simple_filters::DeltaFilterOP > delta_filters() const;
  void delta_filters( utility::vector1< protocols::simple_filters::DeltaFilterOP > const d );
  bool rtmin() const;
  void rtmin( bool const b );
  bool shuffle_order() const;
  void shuffle_order( bool const b );
  bool diversify() const;
  void diversify( bool const b );
  bool incl_nonopt() const;
  void incl_nonopt( bool const b );

private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	protocols::moves::MoverOP relax_mover_;
	utility::vector1< std::string > sample_types_;
	bool dump_pdb_;
	bool dump_table_;
	bool parallel_;
	protocols::filters::FilterOP stopping_condition_; // dflt NULL ; if defined, stops greedy optimization when the filter's apply evaluates to true;
	utility::vector1< std::pair< core::Size, utility::vector1<
			std::pair< core::chemical::AA, utility::vector1< core::Real > > > > > seqpos_aa_vals_vec_;
	utility::vector1< core::Real > filter_deltas_;
	core::pose::Pose ref_pose_;
	utility::vector1< core::pose::Pose > pfront_poses_;
	utility::vector1< utility::vector1< core::Real > > pfront_poses_filter_vals_;
	utility::vector1< utility::vector1< core::Size > > pfront_poses_filter_ranks_;
	core::Size nstruct_iter_;
	bool stop_before_condition_;
	bool skip_best_check_;
	utility::vector1<protocols::simple_filters::DeltaFilterOP> reset_delta_filters_;
	bool rtmin_; //dflt false; should we rtmin after packing?
	core::Real design_shell_;//dflt -1 to only allow pointmutations, higher allows suroundings to be designed as well
	core::Real repack_shell_;
	bool shuffle_order_; //randomize the order that mutations are attempted?
	bool diversify_; //randomize the order that mutations are attempted?
	bool incl_nonopt_; //randomize the order that mutations are attempted?
};


} // design_opt
} // protocols


#endif /*INCLUDED_protocols_design_opt_GreedyOptMutationMover_HH*/
