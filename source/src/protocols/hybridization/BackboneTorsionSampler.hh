// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   BackboneTorsionSampler.hh
///
/// @brief
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_BackboneTorsionSampler_hh
#define INCLUDED_protocols_hybridization_BackboneTorsionSampler_hh

// Unit Headers
//#include <protocols/hybridization/BackboneTorsionSampler.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <protocols/toolbox/task_operations/InterfaceTaskOperation.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.hh>

#include <utility/tag/Tag.fwd.hh>


// Utility Headers
#include <utility/vector1.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers

namespace protocols {
namespace hybridization {

class BackboneTorsionSampler: public protocols::moves::Mover
{
public:
	/// @brief
	///  empty constructor fills values with the expected defaults
	BackboneTorsionSampler();

	//destructor
	~BackboneTorsionSampler() override;

	void init();

	inline void set_scorefunction(core::scoring::ScoreFunctionOP const scorefxn) { scorefxn_ = scorefxn; }

	void local_perturb(core::pose::Pose pose, core::Real max_delta_torsion);

	void perturb(core::pose::Pose & pose,
		core::Size level,
		core::Real max_delta_torsion,
		core::Size local,
		bool rama_biased,
		bool repack,
		bool minimize);

	void apply( core::pose::Pose & pose ) override;

	void task_factory( core::pack::task::TaskFactoryCOP tf );

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & datamap,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	) override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;


	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real increase_cycles_;
	core::Real temperature_;
	bool recover_low_;
	core::Size local_;
	bool dump_snapshots_;
	std::string snapshot_prefix_;
	core::Size snapshot_interval_;

	core::Size n_nested_;
	core::Size perturbed_res_;
	utility::vector1< core::Size > residue_list_;
	protocols::minimization_packing::PackRotamersMoverOP pack_full_repack_;
	core::optimization::AtomTreeMinimizerOP minimizer_;
	//core::optimization::CartesianMinimizerOP minimizer_;

	core::pose::PoseOP native_;

	//AtomTreeMinimizerOP minimizer_;
	core::optimization::MinimizerOptionsOP options_;
	core::kinematics::MoveMap mm_;

	core::pack::task::TaskFactoryCOP task_factory_;
};
} // hybridization
} // protocols

#endif

