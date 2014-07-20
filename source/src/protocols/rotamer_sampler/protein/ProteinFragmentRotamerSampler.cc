// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/protein/ProteinFragmentRotamerSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/protein/ProteinFragmentRotamerSampler.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.protein.ProteinFragmentRotamerSampler" );

///////////////////////////////////////////////////////////////////
//
// Put in for backwards compatibility with
// some 2009 protein stepwise assembly stuff -- not well tested
//        -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	//Constructor
	ProteinFragmentRotamerSampler::ProteinFragmentRotamerSampler( std::string const frag_file,
																									utility::vector1< core::Size > const & slice_res,
																									utility::vector1< core::Size > const & moving_residues	 )
	{
		initialize( frag_file, slice_res, moving_residues );
		set_random( false );
	}

	//Destructor
	ProteinFragmentRotamerSampler::~ProteinFragmentRotamerSampler()
	{}

	/////////////////////////////////////////////////////////////////////////
	void
	ProteinFragmentRotamerSampler::apply( core::pose::Pose & pose, Size const id )
	{
		runtime_assert( id <= size() );
		frame_->apply( id, pose );
	}

	/////////////////////////////////////////////////////////////////////////
	core::Size
	ProteinFragmentRotamerSampler::size() const { return frame_->nr_frags(); }

	/////////////////////////////////////////////////////////////////////////
	void
	ProteinFragmentRotamerSampler::initialize( std::string const frag_file,
																			utility::vector1< core::Size > const & slice_res,
																			utility::vector1< core::Size > const & moving_residues	 ) {

		core::fragment::ConstantLengthFragSetOP fragset =  new core::fragment::ConstantLengthFragSet( 0 /*frag_length ... is reset by reader*/, frag_file );

		if( fragset->max_frag_length() != moving_residues.size() ) {
			utility_exit_with_message( "Number of -moving_res must match frag size!" );
		}
		fragment_set_slice( fragset, slice_res );

		core::fragment::FrameList frames;
		insert_pos_ = moving_residues[ 1 ];
		fragset->frames( insert_pos_, frames );
		frame_ = frames[1];
  }

} //protein
} //rotamer_sampler
} //protocols
