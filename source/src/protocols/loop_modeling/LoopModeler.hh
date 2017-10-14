// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopModeler_HH
#define INCLUDED_protocols_loop_modeling_LoopModeler_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopModeler.fwd.hh>
#include <protocols/loop_modeling/LoopBuilder.fwd.hh>
#include <protocols/loop_modeling/LoopProtocol.fwd.hh>
#include <protocols/loop_modeling/utilities/PrepareForCentroid.fwd.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatom.fwd.hh>
#include <protocols/loop_modeling/LoopModelerTests.fwd.hh> //For friendship

// Core headers
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

// Protocols headers
#include <protocols/filters/Filter.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace loop_modeling {

/// @brief Attempt to find the native structure for one or more loops.
/// @details The typical loop modeling simulation in rosetta has three steps:
/// loop building, centroid refinement, and fullatom refinement.  LoopModeler
/// carries out all three of these steps and allows each to to be enabled,
/// disabled, and otherwise configured.  By default, nothing needs to be
/// specified and a standard loop modeling simulation will be performed.
///
/// Note that this class is a fairly thin wrapper around other LoopMovers.
/// LoopBuilder and LoopProtocol in particular do all the heavy lifting.  The
/// main role of this class is actually to provide a reasonable set of default
/// values, some nice tracer output, and a sophisticated parse_my_tag() method
/// for use with rosetta scripts.
///
/// Note that LoopModeler doesn't implement a proper copy constructor.  (In
/// fact, no LoopMover does.)  This means that if a simulation breaks and
/// nstruct > 1, the remaining simulations will probably break for weird
/// reasons.
class LoopModeler : public LoopMover {

public:

	/// @brief Default constructor.
	LoopModeler();

	/// @brief Default destructor.
	~LoopModeler() override;

	/// @copydoc LoopMover::get_name

	/// @brief Return a shallow copy of this object.
	moves::MoverOP clone() const override;

	/// @copydoc LoopMover::parse_my_tag
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose) override;

protected:

	/// @brief Search for the optimal conformation of the given loops.
	bool do_apply(Pose & pose) override;

public:

	/// @brief Setup the LoopModeler using a minimal configuration.
	void setup_empty_config();

	/// @brief Setup the LoopModeler using the "next-gen KIC" configuration.
	void setup_kic_config();

	/// @brief Setup the LoopModeler using the "KIC with fragments"
	/// configuration.
	void setup_kic_with_fragments_config();

	/// @brief Setup the LoopModeler using the LoopHashKIC configuration
	void setup_loophash_kic_config(bool perturb_sequence);

	/// @brief Return a pointer to the build stage mover.
	LoopBuilderOP build_stage();

	/// @brief Return a pointer to the centroid stage mover.
	LoopProtocolOP centroid_stage();

	/// @brief Return a pointer to the fullatom stage mover.
	LoopProtocolOP fullatom_stage();

	/// @brief Enable the build stage.
	void enable_build_stage();

	/// @brief Disable the build stage.
	void disable_build_stage();

	/// @brief Enable the centroid stage.
	void enable_centroid_stage();

	/// @brief Disable the centroid stage.
	void disable_centroid_stage();

	/// @brief Enable the fullatom stage.
	void enable_fullatom_stage();

	/// @brief Disable the fullatom stage.
	void disable_fullatom_stage();

	/// @brief Get the task factory to be used on the next call to apply().
	/// @details If no task factory has been set, this will raise an exception.
	core::pack::task::TaskFactoryOP get_task_factory();

	/// @brief Get the task factory to be used on the next call to apply().
	/// @details If no task factory has been set, the fallback will be returned.
	core::pack::task::TaskFactoryOP get_task_factory(
		core::pack::task::TaskFactoryOP fallback);

	/// @brief Set the task factory to be used on the next call to apply().
	void set_task_factory(core::pack::task::TaskFactoryOP task_factory);

	/// @brief Set the score function to be used for the fullatom stage.
	void set_fa_scorefxn(core::scoring::ScoreFunctionOP scorefxn);

	/// @brief Set the score function to be used for the centroid stage.
	void set_cen_scorefxn(core::scoring::ScoreFunctionOP scorefxn);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	friend class ::LoopModelerTests;

private:

	LoopBuilderOP build_stage_;
	LoopProtocolOP centroid_stage_;
	LoopProtocolOP fullatom_stage_;
	utilities::PrepareForCentroidOP prepare_for_centroid_;
	utilities::PrepareForFullatomOP prepare_for_fullatom_;

	bool is_build_stage_enabled_;
	bool is_centroid_stage_enabled_;
	bool is_fullatom_stage_enabled_;

};

}
}

#endif
