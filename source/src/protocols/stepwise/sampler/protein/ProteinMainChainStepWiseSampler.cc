// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/protein/ProteinMainChainStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/protein/ProteinMainChainStepWiseSampler.hh>
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sampler.protein.ProteinMainChainStepWiseSampler" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {

//////////////////////////////////////////////////////////////////////////
//constructor!
ProteinMainChainStepWiseSampler::ProteinMainChainStepWiseSampler(
	utility::vector1< core::id::TorsionID > const & which_torsions,
	utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists,
	bool const choose_random /* = false */ ):
	StepWiseSamplerSized(),
	which_torsions_( which_torsions ),
	main_chain_torsion_set_lists_( main_chain_torsion_set_lists )
{
	set_random( choose_random );
}

//////////////////////////////////////////////////////////////////////////
ProteinMainChainStepWiseSampler::ProteinMainChainStepWiseSampler()
{}

//////////////////////////////////////////////////////////////////////////
ProteinMainChainStepWiseSampler::~ProteinMainChainStepWiseSampler()
{}

//////////////////////////////////////////////////////////////////////////
void
ProteinMainChainStepWiseSampler::apply( core::pose::Pose & pose, Size const id )
{
	if ( id > size() ) utility_exit_with_message( "Asked ProteinMainChainStepWiseSampler for another sample but it does not have one!" );
	utility::vector1< Real > const & main_chain_torsion_set_list( main_chain_torsion_set_lists_[ id ] );
	for ( Size i = 1; i <= which_torsions_.size(); i++ ) {
		//std::cout << "SETTING TORSION " << which_torsions_[ i ] << "  to " << main_chain_torsion_set_list[ i ] << std::endl;
		pose.set_torsion( which_torsions_[ i ], main_chain_torsion_set_list[ i ] );
	}
}

///////////////////////////////////////////////////////////////////////////
using namespace core::id;
toolbox::SamplerPlusPlusOP
ProteinMainChainStepWiseSampler::find( TorsionID const & torsion_id ) {
	for ( auto which_torsion_id : which_torsions_ ) {
		if ( which_torsion_id == torsion_id ) return std::dynamic_pointer_cast< ProteinMainChainStepWiseSampler >( shared_from_this() );
	}
	return 0;
}
///////////////////////////////////////////////////////////////////////////


} //protein
} //sampler
} //stepwise
} //protocols
