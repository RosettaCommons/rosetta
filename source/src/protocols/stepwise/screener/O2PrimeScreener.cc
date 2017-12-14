// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/O2PrimeScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/O2PrimeScreener.hh>
#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>
#include <protocols/moves/CompositionMover.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.O2PrimeScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
O2PrimeScreener::O2PrimeScreener( modeler::rna::o2prime::O2PrimePackerOP o2prime_packer ):
	SampleApplier( o2prime_packer->pose() ),
	o2prime_packer_( o2prime_packer )
{}

//Destructor
O2PrimeScreener::~O2PrimeScreener() = default;

/////////////////////////////////////////
bool
O2PrimeScreener::check_screen(){
	o2prime_packer_->sample_o2prime_hydrogen();
	return true;
}

/////////////////////////////////////////
void
O2PrimeScreener::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
	update_mover->add_mover( o2prime_packer_ ); // will copy over 2'-OH solutions
	restore_mover->add_mover( nullptr ); // does nothing to restore -- this is the original choice in SWA
}

} //screener
} //stepwise
} //protocols
