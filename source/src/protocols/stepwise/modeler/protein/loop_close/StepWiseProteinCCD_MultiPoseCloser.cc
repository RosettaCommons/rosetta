// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_MultiPoseCloser.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_MultiPoseCloser.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <core/pose/Pose.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.protein.loop_close.StepWiseProteinCCD_MultiPoseCloser" );

////////////////////////////////////////////////////////////////////////////////////
//
// Wrapper around StepWiseProteinCCD_Closer.cc
//
// A 'temporary' class -- will be deprecated when CCD closer will become
//  a StepWiseScreener (included inside the modeler loop) rather than a sampler
//  in its own right.
//

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace loop_close {

	//Constructor
	StepWiseProteinCCD_MultiPoseCloser::StepWiseProteinCCD_MultiPoseCloser( working_parameters::StepWiseWorkingParametersCOP working_parameters,
																																					sampler::StepWiseSamplerSizedOP sampler ):
		ccd_closer_( new StepWiseProteinCCD_Closer( working_parameters ) ),
		sampler_( sampler ),
		choose_random_( false ),
		num_random_samples_( 20 ),
		max_ntries_( 500 ), // for random modeler.
		ccd_close_res_( 0 )
	{
	}

	//Destructor
	StepWiseProteinCCD_MultiPoseCloser::~StepWiseProteinCCD_MultiPoseCloser()
	{}

  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinCCD_MultiPoseCloser::apply( core::pose::Pose & pose )
	{

		clock_t const time_start( clock() );
		ccd_closer_->set_working_moving_res_list( moving_residues_ );
		ccd_closer_->set_ccd_close_res( ccd_close_res_ );
		ccd_closer_->init( pose );
		main_chain_torsion_sets_.clear();

		Pose pose_save = pose;
		for ( sampler_->reset(); sampler_->not_end(); ++( *sampler_ ) ) {

			if ( choose_random_ && main_chain_torsion_sets_.size() >= num_random_samples_ ) break;
			if ( choose_random_ && ccd_closer_->ntries() >= max_ntries_ ) break;

			sampler_->apply( pose );

			ccd_closer_->apply( pose );
			if ( ccd_closer_->closed_loop() ) {
				TR.Debug << ccd_closer_->which_torsions() << ": " << ccd_closer_->main_chain_torsion_set() << std::endl;
				main_chain_torsion_sets_.push_back( ccd_closer_->main_chain_torsion_set() );
			}

		}

		Size const num_successes = main_chain_torsion_sets_.size();
		std::cout << "CCD closed: " << num_successes << " out of " << ccd_closer_->ntries() << " tries." << std::endl;

		// Kind of silly -- should have at least one output.
		if ( num_successes == 0  ) main_chain_torsion_sets_.push_back( ccd_closer_->grab_main_chain_torsion_set_list( pose ) );

		std::cout << "Total time in StepWiseProteinCCD_Closer: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

		pose = pose_save;

	}

	utility::vector1< id::TorsionID > const &
	StepWiseProteinCCD_MultiPoseCloser::which_torsions() const {
		return ccd_closer_->which_torsions();
	}

	utility::vector1< utility::vector1< core::Real > > const &
	StepWiseProteinCCD_MultiPoseCloser::main_chain_torsion_sets() const {
		return main_chain_torsion_sets_;
	}

} //loop_close
} //protein
} //modeler
} //stepwise
} //protocols
