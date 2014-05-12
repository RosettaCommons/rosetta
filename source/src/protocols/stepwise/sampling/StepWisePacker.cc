// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/StepWisePacker.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/StepWisePacker.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.StepWisePacker" );

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Unified treatment of protein side-chains, RNA 2'-OH packing/trials for stepwise modeling,
//   including smart visualization if there is not information to place the groups.
//
// Also now including 'packing' of terminal phosphates. This is a separate step in which phosphates
//  go on as backbone variants. I also tried a mode where those phosphates were 'side chains' to
//  enable simultaneous optimization with 2'-OH, but code got really slow due to rotamer explosion.
//
//    -- Rhiju, 2014.
//
// To do (maybe):
//   - test instantiation of RNA sugars (would need to update score term so that a single 2'-OH
//      could actually tip the balance and favor instantiation.
//   - develop code to remove/add peptide termini -- N-acetylation and C-terminal methylamidation;
//      currently those variants are just instantiated before sampling and left on even if they
//      don't contact anything.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace sampling {

	//Constructor
	StepWisePacker::StepWisePacker()
	{}

	//Destructor
	StepWisePacker::~StepWisePacker()
	{}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::add_packer_screeners( utility::vector1< StepWiseScreenerOP > & screeners,
																				pose::PoseOP screening_pose ){
		if ( options_->sampler_perform_phosphate_pack() ) screeners_.push_back( new PhosphateScreener( phosphate_sampler_ ) );

		if ( options_->sampler_try_sugar_instantiation() &&
				 working_parameters_->floating_base() )	screeners_.push_back( new SugarInstantiator( *screening_pose_, working_parameters_->working_moving_res() ) );

		if ( options_->sampler_perform_o2prime_pack() )	screeners_.push_back( new O2PrimeScreener( o2prime_packer_ ) );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::initialize( pose::Pose & pose ){

		// actually does a prepack:
		if ( options_->sampler_perform_phosphate_pack() ){
			phosphate_sampler_ = new phosphate::MultiPhosphateSampler( pose, moving_res_list_ );
			runtime_assert(  moving_partition_res_ == phosphate_sampler_->moving_partition_res()/*determined inside*/ );
		}

		if ( o2prime_legacy_mode_ ) {
			if ( options_->sampler_perform_o2prime_pack() ) {
				remove_virtual_O2Prime_hydrogen( pose ); // i.e., instantiate all 2'-OH.
				o2prime_packer_ = new o2prime::O2PrimePacker( pose, scorefxn_, moving_partition_res_ );
			} else {
				add_virtual_O2Prime_hydrogen( pose );
			}
		}

		initialize_protein_packer( pose );
	}

	///////////////////////////////////////////////////////////////////////////
	void
	StepWisePacker::reset( pose::Pose const & pose ){
		phosphate_sampler_->reset( pose );
	}

} //sampling
} //stepwise
} //protocols
