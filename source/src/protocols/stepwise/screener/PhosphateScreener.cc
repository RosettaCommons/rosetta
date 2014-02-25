// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/PhosphateScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/PhosphateScreener.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/moves/CompositionMover.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.PhosphateScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	PhosphateScreener::PhosphateScreener( sampling::rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler ):
		SampleApplier( phosphate_sampler->pose() ),
		phosphate_sampler_( phosphate_sampler ),
		phosphate_sampler_for_restoration_( phosphate_sampler->clone() ) // will not change, and will allow restoration of phosphate.
	{}

	//Destructor
	PhosphateScreener::~PhosphateScreener()
	{}

	/////////////////////////////////////////
	bool
	PhosphateScreener::check_screen(){
		phosphate_sampler_->sample_phosphates();
		return true;
	}

	/////////////////////////////////////////
	void
	PhosphateScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
		update_mover->add_mover( phosphate_sampler_ );

		restore_mover->add_mover( 0 ); // choice in ClassicResidueSampler -- not sure if we should keep it though.
		//		phosphate_sampler_for_restoration_->set_phosphate_move_list( phosphate_sampler_->phosphate_move_list() );
		//		restore_mover->add_mover( phosphate_sampler_for_restoration_ );
	}

} //screener
} //stepwise
} //protocols
