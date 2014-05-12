// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/legacy/StepWiseLegacyConnectionSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/legacy/StepWiseLegacyConnectionSampler.hh>
#include <protocols/stepwise/sampling/protein/legacy/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/sampling/rna/legacy/connection/StepWiseRNA_ConnectionSampler.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampling/working_parameters/util.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.StepWiseLegacyConnectionSampler" );

using namespace core;
using namespace protocols::stepwise::sampling::rna;
using namespace protocols::stepwise::sampling::protein;

namespace protocols {
namespace stepwise {
namespace sampling {

	//Constructor
	StepWiseLegacyConnectionSampler::StepWiseLegacyConnectionSampler( utility::vector1< pose::PoseOP > & pose_list,
																												working_parameters::StepWiseWorkingParametersCOP working_parameters,
																												StepWiseModelerOptionsCOP modeler_options,
																												core::scoring::ScoreFunctionCOP scorefxn):
		pose_list_( pose_list ), // where work will be saved.
		working_parameters_( working_parameters ),
		modeler_options_( modeler_options ),
		scorefxn_( scorefxn ),
		skip_sampling_( false ),
		working_obligate_pack_res_( working_parameters->working_moving_res_list() )
	{}

	//Destructor
	StepWiseLegacyConnectionSampler::~StepWiseLegacyConnectionSampler()
	{}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::apply( pose::Pose & pose ){
		if ( is_protein( pose, working_parameters_->working_moving_res_list() ) ){
			do_protein_residue_sampling( pose );
		} else {
			do_rna_residue_sampling( pose );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::do_protein_residue_sampling( core::pose::Pose & pose ) {

		StepWiseProteinConnectionSampler stepwise_protein_sampler( working_parameters_ );

		stepwise_protein_sampler.set_options( modeler_options_ );
		stepwise_protein_sampler.set_input_streams( input_streams_ );
		stepwise_protein_sampler.set_scorefxn( pack_scorefxn_ );
		stepwise_protein_sampler.set_skip_sampling( skip_sampling_ ); // this means don't do backbone sampling, but do packing.

		stepwise_protein_sampler.apply( pose );

		pose_list_                 = stepwise_protein_sampler.get_pose_list();
		working_obligate_pack_res_ = stepwise_protein_sampler.working_obligate_pack_res();
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::do_rna_residue_sampling( pose::Pose & pose ) {

		using namespace checker;

		if ( skip_sampling_ || working_parameters_->working_moving_res_list().size() == 0 ) {
			add_to_pose_list( pose_list_, pose, "input_pose" );
			return;
		}

		connection::StepWiseRNA_ConnectionSampler stepwise_rna_sampler( working_parameters_ );
		stepwise_rna_sampler.set_silent_file ( modeler_options_->silent_file() + "_sampling" );
		stepwise_rna_sampler.set_scorefxn ( scorefxn_ );
		stepwise_rna_sampler.set_options( modeler_options_->get_sampler_options() );

		stepwise_rna_sampler.apply( pose );

		pose_list_ = stepwise_rna_sampler.get_pose_list();
	}


} //sampling
} //stepwise
} //protocols
