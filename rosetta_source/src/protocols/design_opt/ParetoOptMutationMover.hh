// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/design_opt/ParetoOptMutationMover.hh
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_design_opt_ParetoOptMutationMover_hh
#define INCLUDED_protocols_design_opt_ParetoOptMutationMover_hh
#include <protocols/design_opt/ParetoOptMutationMover.fwd.hh>
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

#ifdef PYROSETTA
	#include <protocols/filters/Filter.hh>
#endif


namespace protocols {
namespace design_opt{

class ParetoOptMutationMover : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	ParetoOptMutationMover();
	ParetoOptMutationMover(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		utility::vector1< protocols::filters::FilterOP > filters,
		utility::vector1< std::string > sample_types,
		bool dump_pdb = false,
		bool dump_table = false,
		bool parallel = false,
		core::Size diversify_lvl = core::Size( 1 ), 
		protocols::filters::FilterOP stopping_condition = protocols::filters::FilterOP( NULL )
	);

	bool pose_coords_are_same( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );
	void filter_seqpos_pareto_opt_ptmuts();
	void clear_cached_data();
	void calc_pfront_poses_filter_ranks();
	void dump_scoring_table( std::string filename, core::pose::Pose const & ref_pose ) const;
	void apply( Pose & pose );
	protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new ParetoOptMutationMover ); }

	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void add_filter(
		protocols::filters::FilterOP,
		std::string const sample_type
	);
	virtual ~ParetoOptMutationMover();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	utility::vector1< protocols::filters::FilterOP > filters() const;
	void filters( utility::vector1< protocols::filters::FilterOP > filters );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	bool dump_pdb() const;
	bool dump_table() const;
	bool parallel() const;
	void dump_pdb( bool const dump_pdb );
	void dump_table( bool const dump_table );
	void parallel( bool const parallel );
	utility::vector1< std::string > sample_types() const;
	void sample_types( utility::vector1< std::string > const sample_types );
	core::Size diversify_lvl() const;
	void diversify_lvl( core::Size const diversify_lvl );
	void stopping_condition( protocols::filters::FilterOP f ){ stopping_condition_ = f; }
	protocols::filters::FilterOP stopping_condition() const{ return stopping_condition_; }
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	protocols::moves::MoverOP relax_mover_;
	utility::vector1< std::string > sample_types_;
	core::Size diversify_lvl_;
	bool dump_pdb_;
	bool dump_table_;
	bool parallel_;
	protocols::filters::FilterOP stopping_condition_; // dflt NULL ; if defined, stops greedy optimization when the filter's apply evaluates to true;
	utility::vector1< std::pair< core::Size, utility::vector1<
			std::pair< core::chemical::AA, utility::vector1< core::Real > > > > > seqpos_aa_vals_vec_;
	core::pose::Pose ref_pose_;
	utility::vector1< core::pose::Pose > pfront_poses_;
	utility::vector1< utility::vector1< core::Real > > pfront_poses_filter_vals_;
	utility::vector1< utility::vector1< core::Size > > pfront_poses_filter_ranks_;
	core::Size nstruct_iter_;
};


} // design_opt
} // protocols


#endif /*INCLUDED_protocols_design_opt_ParetoOptMutationMover_HH*/
