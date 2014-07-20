// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/JumpRotamerSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rotamer_sampler/JumpRotamerSampler.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.JumpRotamerSampler" );

///////////////////////////////////////////////////////////////////
//
//  Bare bones JumpRotamerSampler -- a more versatile version for general
//    rotation/translation searches, which
//    is also less sensitive to jump atoms, allows access to the 'Stub',
//    etc. is in rigid_body/RigidBodyRotamerSampler.
//
//        -- rhiju, 2014
///////////////////////////////////////////////////////////////////

namespace protocols {
namespace rotamer_sampler {

	//Constructor
	JumpRotamerSampler::JumpRotamerSampler():
		which_jump_( 0 )
	{
		set_random( false );
	}

	//Constructor
	JumpRotamerSampler::JumpRotamerSampler(	Size const which_jump,
														utility::vector1< core::kinematics::Jump > const & jumps,
														bool const choose_random /* = false */ ):
		which_jump_( which_jump ),
		jumps_( jumps )
	{
		set_random( choose_random );
	}

	//Destructor
	JumpRotamerSampler::~JumpRotamerSampler()
	{}

  //////////////////////////////////////////////////////////////////////////
	void
	JumpRotamerSampler::apply( core::pose::Pose & pose, Size const id )
	{
		runtime_assert( id <= size() );
		pose.set_jump( which_jump_, jumps_[ id ] );
	}


} //rotamer_sampler
} //protocols
