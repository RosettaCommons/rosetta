// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/design_opt/PointMutationCalculator.hh
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_design_opt_PointMutationCalculator_hh
#define INCLUDED_protocols_design_opt_PointMutationCalculator_hh
#include <protocols/design_opt/PointMutationCalculator.fwd.hh>
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

#include <utility/vector1.hh>

#ifdef PYROSETTA
	#include <protocols/filters/Filter.hh>
#endif


namespace protocols {
namespace design_opt {

class PointMutationCalculator : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
public:
	PointMutationCalculator();
	PointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		utility::vector1< protocols::filters::FilterOP > filters,
		utility::vector1< std::string > sample_types,
		bool dump_pdb = false,
		bool rtmin = false,
		bool parallel = false,
		core::Real design_shell = -1.0,
		core::Real repack_shell = 8.0
	);
	PointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::moves::MoverOP relax_mover,
		protocols::filters::FilterOP filter,
		std::string sample_type = "low",
		bool dump_pdb = false,
		bool rtmin = false,
		bool parallel = false,
		core::Real design_shell = -1.0,
		core::Real repack_shell = 8.0
	);
	virtual ~PointMutationCalculator();

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
		utility::vector1< core::Real > & vals
	);

	void calc_point_mut_filters( Pose const & start_pose,
		utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, utility::vector1< core::Real > > > > > & seqpos_aa_vals_vec );
	void calc_point_mut_filters( Pose const & start_pose,
		utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, core::Real > > > > & seqpos_aa_val_vec );
	protocols::design_opt::PointMutationCalculatorOP clone() const;

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
	utility::vector1< std::string > sample_types() const;
	void sample_types( utility::vector1< std::string > const & sample_types );
	void rtmin( bool const r );
	bool rtmin() const;
	void parallel( bool const r );
	bool parallel() const;
	void set_design_shell( core::Real dz_shell );
	void set_repack_shell( core::Real rp_shell );

private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	protocols::moves::MoverOP relax_mover_;
	utility::vector1< std::string > sample_types_;
	bool dump_pdb_;
	bool rtmin_; //dflt false; should we rtmin after repack?
	bool parallel_; //parallelize calculator with MPI?
	core::Real design_shell_; // dflt -1 which does not mutate the neighbors
	core::Real repack_shell_;
};


} // design_opt
} // protocols


#endif /*INCLUDED_protocols_design_opt_PointMutationCalculator_HH*/
