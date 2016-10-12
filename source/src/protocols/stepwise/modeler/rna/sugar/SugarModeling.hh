// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/rna/sugar/SugarModeling.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_SugarModeling_HH
#define INCLUDED_protocols_stepwise_rna_SugarModeling_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <core/types.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace sugar {

class SugarModeling: public utility::pointer::ReferenceCount {

public:

	SugarModeling( core::Size const input_moving_res, core::Size const input_reference_res );

	SugarModeling();

	~SugarModeling();

	void
	check_compatibility( core::Size const nres ) const;

	void
	set_base_and_pucker_state( core::pose::Pose const & pose, working_parameters::StepWiseWorkingParametersCOP const & WP );

	SugarModeling &
	operator = ( SugarModeling const & src );

public:

	bool sample_sugar;
	core::Size moving_res;
	core::Size reference_res;
	bool is_prepend;
	core::Size bulge_res;
	core::Size bulge_suite;
	core::Size five_prime_chain_break;
	core::chemical::rna::PuckerState moving_res_pucker_state;
	core::chemical::rna::PuckerState bulge_res_pucker_state;
	core::chemical::rna::ChiState moving_res_base_state;
	core::chemical::rna::ChiState bulge_res_base_state;
	utility::vector1< core::pose::PoseOP > pose_list; //pose_data_list of possible sugar conformations.

};


} //sugar
} //rna
} //modeler
} //stepwise
} //protocols

#endif
