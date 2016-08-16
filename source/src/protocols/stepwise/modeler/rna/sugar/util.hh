// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/sugar/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_VirtualSugarUtil_HH
#define INCLUDED_protocols_stepwise_rna_VirtualSugarUtil_HH

#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.fwd.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <map>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace core::pose;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

void
minimize_all_sampled_floating_bases( core::pose::Pose & viewer_pose,
	utility::vector1< SugarModeling > const & modeling_list,
	utility::vector1< PoseOP > & pose_data_list,
	core::scoring::ScoreFunctionCOP const & modeler_scorefxn,
	working_parameters::StepWiseWorkingParametersCOP const & working_parameters,
	bool const virtual_sugar_is_from_prior_step = true );

bool
is_sugar_virtual( core::pose::Pose const & pose, core::Size const sugar_res, core::Size const bulge_res,
	utility::vector1< Size > & bulge_residues_to_virtualize );

bool
is_sugar_virtual( core::pose::Pose const & pose, core::Size const previous_moving_res, core::Size const previous_bulge_res );

void
copy_bulge_res_and_sugar_torsion( SugarModeling const & sugar_modeling, core::pose::Pose & pose, core::pose::Pose const & template_pose,
	bool instantiate_sugar = false );

void
modeler_starting_pose_data_list( utility::vector1< PoseOP > & starting_pose_data_list,
	utility::vector1< SugarModeling > const & SugarModeling_list,
	core::pose::Pose const & pose );


std::map< Size, Size > const
get_reference_res_for_each_virtual_sugar_without_fold_tree( pose::Pose const & pose, Size const moving_suite /*cannot place jump across partititions*/ );

std::map< Size, Size > const
get_reference_res_for_each_virtual_sugar_based_on_fold_tree( pose::Pose const & pose );

Size
get_reference_res_for_virtual_sugar_based_on_fold_tree( pose::Pose const & pose, Size const n );

utility::vector1< Size >
get_possible_reference_res_list_from_pose_without_fold_tree( Size const virtual_sugar_res,
	pose::Pose const & pose,
	Size const  moving_suite /* for old-school poses without jumps already setup*/);


Size
look_for_jumps( Size const n, pose::Pose const & pose, bool const force_upstream );

Size
look_for_jumps_to_previous( Size const virtual_sugar_res,
	pose::Pose const & pose,
	bool const force_upstream );

Size
look_for_jumps_to_next( Size const virtual_sugar_res,
	pose::Pose const & pose,
	bool const force_upstream );

Size
look_for_non_jump_reference_to_previous( Size const virtual_sugar_res,
	pose::Pose const & pose,
	Size const moving_suite );

Size
look_for_non_jump_reference_to_next( Size const virtual_sugar_res,
	pose::Pose const & pose,
	Size const moving_suite );

sampler::copy_dofs::ResidueAlternativeSetOP convert_sugar_modeling_to_residue_alternative_set( SugarModeling const & sugar_modeling );

VirtualSugarJustInTimeInstantiatorOP
instantiate_any_virtual_sugars( pose::Pose & pose,
	working_parameters::StepWiseWorkingParametersCOP working_parameters,
	core::scoring::ScoreFunctionCOP scorefxn,
	options::StepWiseModelerOptionsCOP options );

utility::vector1< bool >
detect_sugar_contacts( pose::Pose const & pose );

bool
detect_sugar_contacts( pose::Pose const & pose, Size const moving_res,
	Distance const o2prime_contact_distance_cutoff_ = 3.2 /*hydrogen bond*/ );

} //sugar
} //rna
} //modeler
} //stepwise
} //protocols

#endif
