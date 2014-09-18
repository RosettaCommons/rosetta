// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/protein/StepWiseProteinMultiPosePacker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/protein/StepWiseProteinMultiPosePacker.hh>
#include <protocols/stepwise/modeler/protein/StepWisePacker.hh>
#include <protocols/stepwise/modeler/protein/checker/ProteinAtrRepChecker.hh>
#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.legacy.modeler.protein.StepWiseProteinMultiPosePacker" );

using namespace ObjexxFCL; // AUTO USING NS
using namespace core;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {
namespace protein {

	//Constructor
	StepWiseProteinMultiPosePacker::StepWiseProteinMultiPosePacker( StepWisePackerOP packer,
																																	sampler::StepWiseSamplerSizedOP sampler ):
		packer_( packer ),
		sampler_( sampler ),
		prepack_( false ),
		atr_rep_check_( false ),
		choose_random_( 0 ),
		num_random_samples_( 20 ),
		max_ntries_( 500 ) // for random modeler.
	{}

	//Destructor
	StepWiseProteinMultiPosePacker::~StepWiseProteinMultiPosePacker()
	{}

	/////////////////////
	std::string
	StepWiseProteinMultiPosePacker::get_name() const {
		return "StepWiseProteinMultiPosePacker";
	}

  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinMultiPosePacker::apply( core::pose::Pose & pose )
	{
		clock_t const time_start( clock() );

		if ( prepack_ ) packer_->do_prepack( pose );
		// atr_rep_checker must be initialized from pose *after* prepack
		if ( atr_rep_check_ ) atr_rep_checker_ = new checker::ProteinAtrRepChecker( pose, moving_res_list_ );
		packer_->initialize( pose );
		sample_residues( pose );

		std::cout << "Total time in StepWiseProteinMultiPosePacker: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " and found " << pose_list_.size() << " samples." << std::endl;
	}


	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinMultiPosePacker::sample_residues( core::pose::Pose & pose )
	{

		 using namespace core::pose;
		 using namespace protocols::stepwise;

		 Pose pose_screen = pose;
		 Size k( 0 );

		 for ( sampler_->reset(); sampler_->not_end(); ++( *sampler_ ) ) {

			 if ( choose_random_ && pose_list_.size() >= num_random_samples_ ) break;
			 if ( choose_random_  && k >= max_ntries_ ) break;
			 k++;

			 sampler_->apply( pose );

			 if ( pose_screener_  && !pose_screener_->check_screen() ) continue;

			 if ( atr_rep_checker_ ){
				 sampler_->apply( pose_screen );
				 if ( !atr_rep_checker_->check_screen( pose_screen ) ) continue;
			 }

			 packer_->apply( pose );
			 pose_list_.push_back( pose.clone() );

		 }

		 //Nothing found? At least produce one pose...
		 if ( pose_list_.size() == 0 ) {
			 TR << "Warning -- nothing passed filters -- outputting a placeholder pose." << std::endl;;
			 pose_list_.push_back( pose.clone() );
		 }

		 sort( pose_list_.begin(), pose_list_.end(), sort_pose_by_score );

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinMultiPosePacker::set_prepack_based_on_moving_res_list( utility::vector1< Size > const & moving_res_list ){
		prepack_ = true;
		moving_res_list_ = moving_res_list; /* also used to set up atr_rep_check */
		packer_->set_moving_res_list( moving_res_list );
	}



} //protein
} //modeler
} //legacy
} //stepwise
} //protocols
