// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Frank DiMaio
#include <protocols/rbsegment_relax/FragInsertAndAlignMover.hh>

// Package headers
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/conformation/util.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>


#include <protocols/frags/RMSVallData.hh>

// C++ Headers
#include <map>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

#include <utility/pointer/ReferenceCount.hh>


#include <string>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace rbsegment_relax {

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "FragInsertAndAlignMover" );


FragInsertAndAlignMover::FragInsertAndAlignMover() {}
FragInsertAndAlignMover::~FragInsertAndAlignMover() {}


FragInsertAndAlignMover::FragInsertAndAlignMover(
	utility::vector1< RBSegment > const &rbsegs,
	core::pose::Pose const &ref_pose,
	core::Real randomness/*=0.0*/ ) {
	initialize_rb_fragments(rbsegs,ref_pose,randomness );
}


/// @brief Initialize fragment library
/// @brief     Loads RMS vall data, sets up frame set
void FragInsertAndAlignMover::initialize_rb_fragments(
	utility::vector1< RBSegment > const &rbsegs,
	core::pose::Pose const &ref_pose,
	core::Real randomness/*=0.0*/ )
{
	using namespace basic::options;

	if ( !option[ OptionKeys::loops::vall_file ].user() ) {
		TR.Error << "FragInsertAndAlignMover::initialize_rb_fragments()  with no vall!  Specify location with -loops::vall" << std::endl;
		return;
	}

	if ( rbsegs.size() == 0 ) {
		TR.Error << "FragInsertAndAlignMover::initialize_rb_fragments() called with empty RBSeglist! continuing ..." << std::endl;
		return;
	}

	// declare as static so we only load once
	static protocols::frags::RMSVallData rms_vall( option[ OptionKeys::loops::vall_file ] );

	std::string input_seq = ref_pose.sequence();
	//core::Size res_ctr = 1;

	frames_.clear();

	// TO DO: COMPOUND SEGMENT SUPPORT(??)
	for ( int i =  1; i <= (int)rbsegs.size(); ++i ) {
		if ( rbsegs[i].isCompound() ) {
			TR.Warning << "[ WARNING ]  FragInsertAndAlignMover::initialize_rb_fragments() undefined for compound segments! continuing..." << std::endl;
			continue;
		}

		int frag_start = rbsegs[i][1].start();
		int frag_size  = rbsegs[i][1].end() - frag_start + 1;

		// add a new frame
		frames_.push_back( core::fragment::FrameOP( new core::fragment::Frame( frag_start, frag_size ) ) );

		utility::vector1< numeric::xyzVector< core::Real> > cas( frag_size );
		for ( int k=0; k<frag_size; ++k ) {
			cas[k+1] = ref_pose.residue(frag_start+k).atom("CA").xyz();
		}

		std::string frag_seq = input_seq.substr( rbsegs[i][1].start()-1, frag_size );
		rms_vall.get_frags(
			200, cas, frag_seq,  rbsegs[i][1].char_type(), frames_[frames_.size()], randomness
		);
	}
	TR.Debug << "Read " << frames_.size() << " frames" << std::endl;
}


void FragInsertAndAlignMover::apply( core::pose::Pose & pose ) {
	if ( frames_.size() == 0 ) {
		TR.Warning << "[ WARNING ]  FragInsertAndAlignMover::apply()called with no frames defined! continuing..." << std::endl;
		return;
	}

	// pick a random rbseg
	int idx = numeric::random::random_range( 1, frames_.size() );
	apply( pose,idx );
}

std::string
FragInsertAndAlignMover::get_name() const {
	return "FragInsertAndAlignMover";
}


void FragInsertAndAlignMover::apply( core::pose::Pose & pose, int idx, bool idealize ) {
	core::Size start = frames_[idx]->start(),len = frames_[idx]->length();

	TR.Debug << "FragInsertAndAlignMover::apply() called [" << start << " , " << start+len-1 << "]" << std::endl;

	// grab coords
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, len );
	numeric::xyzVector< core::Real > com1(0,0,0);
	for ( int i=0; i<(int)len; ++i ) {
		numeric::xyzVector< core::Real > x_i = pose.residue(start+i).atom("CA").xyz();  // CA
		com1 += x_i;
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i+1) = x_i[j];
		}
	}
	com1 /= len;
	for ( int i=0; i<(int)len; ++i ) {
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i+1) -= com1[j];
		}
	}

	//core::pose::Pose pose_copy = pose;

	if ( idealize ) {
		for ( int i=0; i<(int)len; ++i ) {
			core::conformation::idealize_position(start+i, pose.conformation());
		}
	}

	// insert frag
	core::Size toget = numeric::random::random_range( 1, frames_[idx]->nr_frags() );
	frames_[idx]->apply( toget, pose );

	// grab new coords
	ObjexxFCL::FArray2D< core::Real > final_coords( 3, len );
	numeric::xyzVector< core::Real > com2(0,0,0);
	for ( int i=0; i<(int)len; ++i ) {
		numeric::xyzVector< core::Real > x_i = pose.residue(start+i).atom(" CA ").xyz();  // CA
		com2 += x_i;
		for ( int j=0; j<3; ++j ) {
			final_coords(j+1,i+1) = x_i[j];
		}
	}
	com2 /= len;
	for ( int i=0; i<(int)len; ++i ) {
		for ( int j=0; j<3; ++j ) {
			final_coords(j+1,i+1) -= com2[j];
		}
	}

	// get optimal superposition
	// rotate >final< to >init<
	ObjexxFCL::FArray1D< numeric::Real > ww( len, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	numeric::model_quality::findUU( final_coords, init_coords, ww, len, uu, ctx );
	numeric::xyzMatrix< core::Real > R;
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	// apply rotation to ALL atoms
	// x_i' <- = R*x_i + com1;


	for ( Size i = 0; i < len; ++i ) {
		//std::cout << i << " " << j << "  " << pose.xyz(id) << "   " << R * ( pose.xyz(id) - com2) + com1
		for ( Size j = 1; j <= pose.residue_type(start+i).natoms(); ++j ) {
			id::AtomID id( j, start+i );
			pose.set_xyz( id, R * ( pose.xyz(id) - com2) + com1 );
		}
	}
}

/// @brief take a CA-only pose and insert idealized fragments close to the trace
void FragInsertAndAlignMover::bootstrapCATrace( core::pose::Pose & start_pose ) {
	for ( int i=1; i<=(int)frames_.size(); ++i ) {
		apply( start_pose , i , true);
	}
}


}
}
