// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/protein/ProteinFragmentStepWiseSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/protein/ProteinFragmentStepWiseSampler.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.sampler.protein.ProteinFragmentStepWiseSampler" );

///////////////////////////////////////////////////////////////////
//
// Put in for backwards compatibility with
// some 2009 protein stepwise assembly stuff -- not well tested
//        -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace sampler {
namespace protein {

//Constructor
ProteinFragmentStepWiseSampler::ProteinFragmentStepWiseSampler( std::string const & frag_file,
	utility::vector1< core::Size > const & slice_res,
	utility::vector1< core::Size > const & moving_residues  )
{
	initialize( frag_file, slice_res, moving_residues );
	set_random( false );
}

//Destructor
ProteinFragmentStepWiseSampler::~ProteinFragmentStepWiseSampler() = default;

/////////////////////////////////////////////////////////////////////////
void
ProteinFragmentStepWiseSampler::apply( core::pose::Pose & pose, Size const id )
{
	runtime_assert( id <= size() );
	frame_->apply( id, pose );
}

/////////////////////////////////////////////////////////////////////////
core::Size
ProteinFragmentStepWiseSampler::size() const { return frame_->nr_frags(); }

/////////////////////////////////////////////////////////////////////////
void
ProteinFragmentStepWiseSampler::initialize( std::string const & frag_file,
	utility::vector1< core::Size > const & slice_res,
	utility::vector1< core::Size > const & moving_residues  ) {

	core::fragment::ConstantLengthFragSetOP fragset( new core::fragment::ConstantLengthFragSet( 0 /*frag_length ... is reset by reader*/, frag_file ) );

	if ( fragset->max_frag_length() != moving_residues.size() ) {
		utility_exit_with_message( "Number of -moving_res must match frag size!" );
	}
	fragment_set_slice( fragset, slice_res );

	core::fragment::FrameList frames;
	insert_pos_ = moving_residues[ 1 ];
	fragset->frames( insert_pos_, frames );
	frame_ = frames[1];
}

} //protein
} //sampler
} //stepwise
} //protocols
