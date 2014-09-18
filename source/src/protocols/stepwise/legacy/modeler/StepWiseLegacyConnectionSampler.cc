// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/modeler/StepWiseLegacyConnectionSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/legacy/modeler/StepWiseLegacyConnectionSampler.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinConnectionSampler.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/legacy/modeler/rna/connection/StepWiseRNA_ConnectionSampler.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/working_parameters/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.modeler.StepWiseLegacyConnectionSampler" );

using namespace core;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler::protein;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {

	//Constructor
	StepWiseLegacyConnectionSampler::StepWiseLegacyConnectionSampler( utility::vector1< pose::PoseOP > & pose_list,
																												stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
																												options::StepWiseModelerOptionsCOP options,
																												core::scoring::ScoreFunctionCOP scorefxn):
		pose_list_( pose_list ), // where work will be saved.
		working_parameters_( working_parameters ),
		options_( options ),
		scorefxn_( scorefxn ),
		skip_modeler_( false ),
		working_obligate_pack_res_( working_parameters->working_moving_res_list() )
	{}

	//Destructor
	StepWiseLegacyConnectionSampler::~StepWiseLegacyConnectionSampler()
	{}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::apply( pose::Pose & pose ){
		if ( is_protein( pose, working_parameters_->working_moving_res_list() ) ){
			do_protein_residue_modeler( pose );
		} else {
			do_rna_residue_modeler( pose );
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::do_protein_residue_modeler( core::pose::Pose & pose ) {

		StepWiseProteinConnectionSampler stepwise_protein_sampler( working_parameters_ );

		stepwise_protein_sampler.set_options( options_ );
		stepwise_protein_sampler.set_input_streams( input_streams_ );
		stepwise_protein_sampler.set_scorefxn( pack_scorefxn_ );
		stepwise_protein_sampler.set_skip_modeler( skip_modeler_ ); // this means don't do backbone modeler, but do packing.

		stepwise_protein_sampler.apply( pose );

		pose_list_                 = stepwise_protein_sampler.get_pose_list();
		working_obligate_pack_res_ = stepwise_protein_sampler.working_obligate_pack_res();
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	StepWiseLegacyConnectionSampler::do_rna_residue_modeler( pose::Pose & pose ) {

		using namespace checker;

		if ( skip_modeler_ || working_parameters_->working_moving_res_list().size() == 0 ) {
			add_to_pose_list( pose_list_, pose, "input_pose" );
			return;
		}

		connection::StepWiseRNA_ConnectionSampler stepwise_rna_sampler( working_parameters_ );
		stepwise_rna_sampler.set_silent_file ( options_->silent_file() + "_modeler" );
		stepwise_rna_sampler.set_scorefxn ( scorefxn_ );
		stepwise_rna_sampler.set_options( options_->get_sampler_options() );

		stepwise_rna_sampler.apply( pose );

		pose_list_ = stepwise_rna_sampler.get_pose_list();
	}


} //modeler
} //legacy
} //stepwise
} //protocols
