// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteMover.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_monte_carlo_DeleteMover_hh
#define INCLUDED_protocols_stepwise_monte_carlo_DeleteMover_hh

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMoverCreator.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.fwd.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class DeleteMover: public protocols::moves::Mover {
public:

	//constructor
	DeleteMover();

	//destructor
	~DeleteMover();

public:

	using protocols::moves::Mover::apply;

	void
	apply( core::pose::Pose & pose, core::Size const res_to_delete_in_full_model_numbering );

	void
	apply( core::pose::Pose & pose, utility::vector1< core::Size > const & residues_to_delete_in_full_model_numbering );

	/// @brief Apply the minimizer to one pose
	virtual void apply( core::pose::Pose & pose_to_visualize ) override;
	protocols::moves::MoverOP fresh_instance() const override { return DeleteMoverOP( new DeleteMover ); }
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


	bool
	decide_to_keep_pose( core::pose::Pose const & pose ) const;

	void
	remove_singletons_and_update_pose_focus( core::pose::Pose & pose,
		core::pose::PoseOP sliced_out_pose_op,
		bool & keep_remainder_pose,
		bool & keep_sliced_out_pose ) const;

	void
	wipe_out_moving_residues( core::pose::Pose & pose );

	void minimize_after_delete( core::pose::Pose & pose ) const;
	void set_minimize_after_delete( bool const setting ){ minimize_after_delete_ = setting; }

	void set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler );

	void
	set_options( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options );

private:

	protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler_;
	protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options_;
	bool minimize_after_delete_;
	utility::vector1< core::Size > interface_res_;
	// for RosettaScripts style use
	utility::vector1< core::Size > residues_to_delete_in_full_model_numbering_;

};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
