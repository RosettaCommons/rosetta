// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/matdes/MatDesPointMutationCalculator.hh
/// @brief this is a modified version of chris king's PointMutationCalculator with additional functionality that is currently not compatible with all of the ParetoOpt functionality. please note that this has been checked into master in its current state in response to requests from others to use this modified version of chris king's GreedyOptMutationMover. although this is still a somewhat developmental piece of code, it has currently been left in src/protocols/matdes/ to avoid issues with intra-library level dependencies.
/// @author jacob bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_MatDesPointMutationCalculator_hh
#define INCLUDED_protocols_matdes_MatDesPointMutationCalculator_hh
#include <protocols/matdes/MatDesPointMutationCalculator.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/simple_filters/TaskAwareScoreTypeFilter.fwd.hh>

// Stored tasks
#include <protocols/toolbox/task_operations/STMStoredTask.hh>

#include <utility/vector1.hh>

#ifdef PYROSETTA
	#include <protocols/filters/Filter.hh>
#endif


namespace protocols {
namespace matdes {

class MatDesPointMutationCalculator : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
public:
	MatDesPointMutationCalculator();
	MatDesPointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		utility::vector1< protocols::filters::FilterOP > filters,
		utility::vector1< std::pair< std::string, core::Real > > filter_thresholds,
		utility::vector1< protocols::simple_filters::TaskAwareScoreTypeFilterOP > set_task_for_filters,
		utility::vector1< std::string > sample_types,
		bool dump_pdb = false,
		bool rtmin = false,
		bool parallel = false,
		core::Real design_shell = -1.0,
		core::Real repack_shell = 8.0,
		bool force_natro = false
	);
	MatDesPointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		protocols::filters::FilterOP filter,
		std::pair< std::string, core::Real > filter_threshold,
		protocols::simple_filters::TaskAwareScoreTypeFilterOP set_task_for_filter,
		std::string sample_type = "low",
		bool dump_pdb = false,
		bool rtmin = false,
		bool parallel = false,
		core::Real design_shell = -1.0,
		core::Real repack_shell = 8.0,
		bool force_natro = false
	);
	virtual ~MatDesPointMutationCalculator();

	void mutate_and_relax(
		core::pose::Pose & pose,
		core::Size const & resi,
		core::chemical::AA const & target_aa
	);

	void mutate_and_relax(
		core::pose::Pose & pose,
		core::Size const & resi,
		core::chemical::AA const & target_aa,
		protocols::simple_moves::GreenPackerOP green_packer
	);

	void eval_filters(
		core::pose::Pose & pose,
		bool & filter_pass,
		utility::vector1< core::Real > & vals,
		bool always_eval
	);

	void calc_point_mut_filters( Pose const & start_pose,
		utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, utility::vector1< core::Real > > > > > & seqpos_aa_vals_vec );
	void calc_point_mut_filters( Pose const & start_pose,
		utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, core::Real > > > > & seqpos_aa_val_vec );
	protocols::matdes::MatDesPointMutationCalculatorOP clone() const;

	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	utility::vector1< protocols::filters::FilterOP > filters() const;
	void filters( utility::vector1< protocols::filters::FilterOP > filters );
	utility::vector1< std::pair< std::string, core::Real > > filter_thresholds() const;
	void filter_thresholds( utility::vector1< std::pair< std::string, core::Real > > filter_thresholds );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	bool dump_pdb() const;
	void dump_pdb( bool const dump_pdb );
	utility::vector1< std::string > sample_types() const;
	void sample_types( utility::vector1< std::string > const & sample_types );
	void rtmin( bool const r );
	bool rtmin() const;
	void parallel( bool const r );
	bool parallel() const;
	void set_design_shell( core::Real dz_shell );
	void set_repack_shell( core::Real rp_shell );
	utility::vector1< protocols::simple_filters::TaskAwareScoreTypeFilterOP > set_task_for_filters() const;
	void set_task_for_filters( utility::vector1< protocols::simple_filters::TaskAwareScoreTypeFilterOP > set_task_for_filters );

	void force_natro( bool const force_natro );
	bool force_natro() const;
	void stored_task_names( utility::vector1< std::string > stored_task_names );
	utility::vector1< std::string > stored_task_names() const;
	void new_tasks( utility::vector1< core::pack::task::PackerTaskOP > new_tasks);
	utility::vector1< core::pack::task::PackerTaskOP > new_tasks() const;
	core::pose::PoseOP reference_pose() const;
	void reference_pose( core::pose::PoseOP const reference_pose );
	//void set_stored_and_previous_tasks( core::pose::Pose & pose );

private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	utility::vector1< std::pair< std::string, core::Real > > filter_thresholds_;
	protocols::moves::MoverOP relax_mover_;
	utility::vector1< std::string > sample_types_;
	bool dump_pdb_;
	bool rtmin_; //dflt false; should we rtmin after repack?
	bool parallel_; //parallelize calculator with MPI?
	core::Real design_shell_; // dflt -1 which does not mutate the neighbors
	core::Real repack_shell_;
	utility::vector1< protocols::simple_filters::TaskAwareScoreTypeFilterOP > set_task_for_filters_;
	bool force_natro_;
	utility::vector1< std::string > stored_task_names_;
	utility::vector1< core::pack::task::PackerTaskOP > new_tasks_;
	core::pose::PoseOP reference_pose_;
};


} // matdes
} // protocols


#endif /*INCLUDED_protocols_matdes_MatDesPointMutationCalculator_HH*/
