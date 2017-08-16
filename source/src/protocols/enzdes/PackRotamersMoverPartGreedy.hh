// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sagar Khare (khares@uw.edu)

#ifndef INCLUDED_protocols_enzdes_PackRotamersMoverPartGreedy_hh
#define INCLUDED_protocols_enzdes_PackRotamersMoverPartGreedy_hh

// Unit headers
//#include <protocols/simple_moves/PackRotamersMoverPartGreedy.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

/// @brief a mover that packs the side-chains around a given set of target residues in a greedy fashion, and then packs the rest using the sim annealer.

class PackRotamersMoverPartGreedy: public protocols::moves::Mover {

	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	/// @brief default constructors
	PackRotamersMoverPartGreedy(
		ScoreFunctionOP scorefxn,
		PackerTaskOP task,
		utility::vector1 <core::Size> target_residues
	);

	PackRotamersMoverPartGreedy();

	//std::string PackRotamersMoverPartGreedyCreator::mover_name();

	// destructor
	~PackRotamersMoverPartGreedy() override;

	//parser stuff
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &pose ) override;
	void apply( core::pose::Pose &pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	// methods
	void greedy_around(
		core::pose::Pose & pose,
		utility::vector1<core::Size > const & target_res,
		core::pack::task::PackerTaskOP task,
		core::scoring::ScoreFunctionCOP scorefxn
	);

	utility::vector1< core::Size > compute_designable_neighbors(
		core::Size const & position,
		core::pack::task::PackerTaskCOP task,
		core::pose::Pose const & pose
	);

	void update_task_and_neighbors(
		core::Size const & best_neigh,
		core::pack::task::PackerTaskOP task,
		utility::vector1< core::Size > & current_neighbors
	);

	// getters and setters

	void task_factory( core::pack::task::TaskFactoryOP p );
	void task( core::pack::task::PackerTaskOP task );
	void target_residues (utility::vector1< core::Size > & trg_res);

	//choose n best residues interacting with ligand
	utility::vector1<core::Size> choose_n_best(
		core::pose::Pose const & pose,
		core::Size const & n_best
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::scoring::ScoreFunctionOP scorefxn_repack_;
	core::scoring::ScoreFunctionOP scorefxn_repack_greedy_;
	core::scoring::ScoreFunctionOP scorefxn_minimize_;
	core::pack::task::PackerTaskOP task_;
	core::pack::task::TaskFactoryOP task_factory_;
	utility::vector1< core::Size > target_residues_;
	utility::vector1< core::Size > restrict_to_repacking_;
	bool use_cstids_;
	core::Real threshold_;
	std::string cstid_list_;
	core::Size n_best_;
};


} // moves
} // protocols

#endif
