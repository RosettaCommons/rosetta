// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/GraftedStemOptimizer.hh
/// @brief Optimize the CDR Grafted stem
/// @details
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody_GraftedStemOptimizer_hh
#define INCLUDED_protocols_antibody_GraftedStemOptimizer_hh

#include <protocols/antibody/GraftedStemOptimizer.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>


namespace protocols {
namespace antibody {

/// @brief Grafts only one CDR onto a framework
class GraftedStemOptimizer : public protocols::moves::Mover {
public:


	// constructor with arguments
	GraftedStemOptimizer(CDRNameEnum const & cdr_name,
		AntibodyInfoOP antibody_info);


	~GraftedStemOptimizer();

	virtual void
	apply( core::pose::Pose & pose );

	virtual std::string
	get_name() const;


	/// @brief enable benchmark mode
	inline void
	enable_benchmark_mode( bool const & setting ) {
		benchmark_ = setting;
	}

	/// @brief users can pass their own scorefunction, foldtree, movemap, and taskfactory
	void
	set_scorefxn(core::scoring::ScoreFunctionOP setting);

	/// @brief deep optimization will do small_share_ccd
	void set_deep_optimization(bool const & setting) {
		deep_optimization_ = setting;
	}

	core::kinematics::FoldTreeOP
	get_N_C_stems_foldtree(core::pose::Pose const & pose) const;

	core::kinematics::FoldTreeOP
	get_Nstem_foldtree(core::pose::Pose const & pose) const;

	core::kinematics::FoldTreeOP
	get_Cstem_foldtree(core::pose::Pose const & pose) const;

	core::kinematics::MoveMapOP
	get_stem_movemap(core::pose::Pose const & pose,
		std::string const & type,
		bool const & include_nb_sc = false) const;

	core::pack::task::TaskFactoryOP
	get_stem_taskfactory(core::pose::Pose & pose,
		std::string const & type,
		bool const & include_nb_sc = false) const;


	/// @brief stem that was replaced by the extra reesidues at the end of the
	///       loop terminus. Default is "2"
	void
	set_stem_size(core::Size const & setting);

	/// @brief copy ctor
	GraftedStemOptimizer( GraftedStemOptimizer const & rhs );

	/// @brief assignment operator
	GraftedStemOptimizer & operator=( GraftedStemOptimizer const & rhs );


private:
	void
	init();

	void
	setup_protocol(core::pose::Pose & pose);

	void
	initForEqualOperatorAndCopyConstructor(GraftedStemOptimizer & lhs, GraftedStemOptimizer const & rhs);


private:
	mutable Size stem_size_;
	CDRNameEnum cdr_name_;
	AntibodyInfoOP ab_info_;
	bool benchmark_;
	loops::LoopOP cdr_loop_;
	core::scoring::ScoreFunctionOP scorefxn_;
	moves::MonteCarloOP mc_;
	moves::SequenceMoverOP optimize_stems_;
	bool deep_optimization_;


}; // class GraftedStemOptimizer


} // antibody
} // protocols


#endif
