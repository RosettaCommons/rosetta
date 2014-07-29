// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_protein_StepWiseProteinWorkingParametersUtil_HH
#define INCLUDED_protocols_stepwise_modeler_protein_StepWiseProteinWorkingParametersUtil_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseProteinWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace working_parameters {

	StepWiseProteinWorkingParametersOP
	setup_protein_working_parameters_for_swa( utility::vector1< Size > const & moving_res_list,
																						pose::Pose const & pose,
																						pose::PoseCOP native_pose,
																						utility::vector1< Size > const & bridge_res,
																						utility::vector1< Size > const & working_minimize_res );

} //working_parameters
} //protein
} //modeler
} //stepwise
} //protocols

#endif
