// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/FromScratchMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_FromScratchMover_HH
#define INCLUDED_protocols_stepwise_monte_carlo_FromScratchMover_HH

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMoverCreator.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

class FromScratchMover: public protocols::moves::Mover {

public:

	//constructor
	FromScratchMover();

	//destructor
	~FromScratchMover();

public:

	using moves::Mover::apply;

	void
	apply( core::pose::Pose & pose,
		utility::vector1<core::Size> const & residues_to_instantiate_in_full_model_numbering ) const;

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize ) override;
	protocols::moves::MoverOP fresh_instance() const override { return FromScratchMoverOP( new FromScratchMover ); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	void set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler );

private:

	void
	update_full_model_info_and_switch_focus_to_new_pose( core::pose::Pose & pose, core::pose::Pose & new_pose, utility::vector1< core::Size > const & resnum ) const;

	void
	sample_by_swa( core::pose::Pose & pose, core::Size const sample_res ) const;

private:

	protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler_;
	utility::vector1< core::Size > residues_to_instantiate_in_full_model_numbering_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
