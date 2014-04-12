// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ProteinPackScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/ProteinPackScreener.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ProteinPackScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	ProteinPackScreener::ProteinPackScreener( pose::Pose & pose,
																						sampling::protein::StepWiseProteinPackerOP stepwise_packer ):
		pose_( pose ),
		stepwise_packer_( stepwise_packer )
	{}

	//Destructor
	ProteinPackScreener::~ProteinPackScreener()
	{}

	bool
	ProteinPackScreener::check_screen() {
		stepwise_packer_->apply( pose_ );
		return true;
	}

} //screener
} //stepwise
} //protocols
