// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/JumpSRFD.cc
/// @brief  A fragment that changes a RigidBody Transform ( Jump )
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/JumpSRFD.hh>

// Package Headers
#include <core/fragment/Frame.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedStubID.hh>

#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

namespace core {
namespace fragment {

using namespace ObjexxFCL::format;

static thread_local basic::Tracer tr( "core.fragment.JumpSRFD" );
using namespace core::id;
//void DownJumpSRFD::get_stubs( pose::Pose const& pose, Size downstream_res_nr, id::StubID &up_stub, id::StubID &down_stub ) const {

// kinematics::Edge jump_edge = pose.fold_tree().get_residue_edge( downstream_res_nr );
// runtime_assert( jump_edge.is_jump() ); // being paranoid
//work out upstream and downstream residue
// Size res1=jump_edge.start();
//Size res2=jump_edge.stop();
//obsolet code:
// runtime_assert( 0 );
//up_stub = StubID( NamedStubID( upstream_stub_atoms_, res1 ), pose );
// down_stub = StubID( NamedStubID( downstream_stub_atoms_, res2 ), pose );


//  // work out the stubID
//  chemical::ResidueType const& rt1 ( pose.residue_type ( res1 ) );
//  chemical::ResidueType const& rt2 ( pose.residue_type ( res2 ) );

//  id::AtomID b1( rt2.atom_index (upstream_stub_atoms_[ 1 ]) , res2 );
//  id::AtomID b2( rt2.atom_index (upstream_stub_atoms_[ 2 ]) , res2 );
//  id::AtomID b3( rt2.atom_index (upstream_stub_atoms_[ 3 ]) , res2 );
//  up_stub = id::StubID( b1, b2, b3 );

//  id::AtomID a1( rt1.atom_index (downstream_stub_atoms_[ 1 ]) , res1 );
//  id::AtomID a2( rt1.atom_index (downstream_stub_atoms_[ 2 ]) , res1 );
//  id::AtomID a3( rt1.atom_index (downstream_stub_atoms_[ 3 ]) , res1 );
//  down_stub = id::StubID( a1, a2, a3 );

//}


/// @brief for DownJumpSRFD this function should never be called, instead use Frame version
/// @return always false
/// @warning will trigger a false runtime assert
bool DownJumpSRFD::apply( kinematics::MoveMap const &, pose::Pose &, Size const ) const {
	runtime_assert( 0 );
	return false;
}


bool DownJumpSRFD::apply( pose::Pose& pose, Size ipos, Frame const& frame ) const {
	//at this point we are the downstream partner
	runtime_assert( ipos > 1 );
	Size const upstream_resnr ( frame.seqpos( ipos-1 ) );
	Size const downstream_resnr( frame.seqpos( ipos ) );
	tr.Trace << "DownJumpSRFD applied for jump " << upstream_resnr << " ---> "<< downstream_resnr << std::endl;

	//check if there is a jump!
	bool has_jump ( pose.fold_tree().is_jump_point( upstream_resnr ) && pose.fold_tree().is_jump_point( downstream_resnr ) );
	if ( !has_jump ) {
		tr.Warning << "WARNING: DownJump-Frag could not be applied to current pose... need to have a jump at "
			<< upstream_resnr << " " << downstream_resnr << std::endl;
		return false;
	}

	// do the dof-change
	id::StubID up_stub( pose::named_stub_id_to_stub_id( NamedStubID( upstream_stub_atoms_,   upstream_resnr   ), pose ) );
	id::StubID down_stub( pose::named_stub_id_to_stub_id( NamedStubID( downstream_stub_atoms_, downstream_resnr ), pose ) );
	tr.Trace << "apply stub transform to " << std::endl;
	tr.Trace << "up stub: " << NamedStubID( upstream_stub_atoms_, upstream_resnr ) << std::endl;
	tr.Trace << "down stub: " << NamedStubID( downstream_stub_atoms_, downstream_resnr ) << std::endl;
	pose.conformation().set_stub_transform( up_stub, down_stub, rt_ );
	return true;
}


bool DownJumpSRFD::apply(
	kinematics::MoveMap const & movemap,
	pose::Pose & pose,
	Size const intra_frame_pos,
	Frame const & frame
) const {

	Size const upstream_resnr( frame.seqpos( intra_frame_pos - 1 ) );
	Size const downstream_resnr( frame.seqpos( intra_frame_pos ) );

	if ( movemap.get_jump( upstream_resnr, downstream_resnr ) ) {
		return apply( pose, intra_frame_pos, frame );
	} else {
		return false;
	}
}


bool DownJumpSRFD::steal( pose::Pose const& pose, Size ipos, Frame const& frame ) {
	//at this point we are the downstream partner
	runtime_assert( ipos > 1 );
	Size const upstream_resnr ( frame.seqpos( ipos-1 ) );
	Size const downstream_resnr( frame.seqpos( ipos ) );


	id::StubID up_stub( pose::named_stub_id_to_stub_id( NamedStubID( upstream_stub_atoms_,   upstream_resnr   ), pose ) );
	id::StubID down_stub( pose::named_stub_id_to_stub_id( NamedStubID( downstream_stub_atoms_, downstream_resnr ), pose ) );
	rt_ = pose.conformation().get_stub_transform( up_stub, down_stub );
	return true; //can something go wrong ? check has_torsion() ?
}

bool DownJumpSRFD::is_compatible( SingleResidueFragData const& aSRFD) const {
	//DownJumpSRFD const* ptr = dynamic_cast< DownJumpSRFD const* > ( & aSRFD );
	//if ( ptr ) {
	//  return true;
	//}
	return dynamic_cast< DownJumpSRFD const* > ( & aSRFD );
}


bool UpJumpSRFD::is_compatible( SingleResidueFragData const& aSRFD) const {
	//UpJumpSRFD const* ptr = dynamic_cast< UpJumpSRFD const* > ( & aSRFD );
	//if ( ptr ) {
	//  return true;
	//}
	return dynamic_cast< UpJumpSRFD const* > ( & aSRFD );
}

bool DownJumpSRFD::is_applicable(
	kinematics::MoveMap const& move_map,
	Size  ipos,
	Frame const& frame
) const {
	runtime_assert( ipos > 1 ); //UpJump DownJump:  must be at least second in Frame
	Size const upstream_resnr ( frame.seqpos( ipos-1 ) );
	Size const downstream_resnr( frame.seqpos( ipos ) );
	return move_map.get_jump( upstream_resnr, downstream_resnr );
}


void DownJumpSRFD::set_stub_atoms( AtomList downstream, AtomList upstream  ) {
	upstream_stub_atoms_ = upstream;
	downstream_stub_atoms_ = downstream;
}

void DownJumpSRFD::set_standard_stub_atoms() {
	downstream_stub_atoms_.reserve(4);
	downstream_stub_atoms_.push_back( "CA" );
	downstream_stub_atoms_.push_back( "N" );
	downstream_stub_atoms_.push_back( "CA" );
	downstream_stub_atoms_.push_back( "C" );
	upstream_stub_atoms_ = downstream_stub_atoms_;
}


void DownJumpSRFD::show( std::ostream &out ) const {
	Parent::show( out );

	//only the downstream_partner has to show his cards
	out << " "<< rt_;
	out << " UPSTUB ";
	for ( AtomList::const_iterator it = upstream_stub_atoms_.begin(),
			eit = upstream_stub_atoms_.end(); it!=eit; ++it ) {
		out << A(3, (*it));
	}
	if ( upstream_stub_atoms_ != downstream_stub_atoms_ ) {
		out << " DOWNSTUB ";
		for ( AtomList::const_iterator it = downstream_stub_atoms_.begin(),
				eit = downstream_stub_atoms_.end(); it!=eit; ++it ) {
			out << A(3, (*it));
		}
	}
}

void DownJumpSRFD::read( std::istream &in ) {
	Parent::read_data( in );
	in >> rt_;
	tr.Trace << "read RT " << rt_ << " status of stream: " << (in.fail() ? "FAIL" : "GOOD" ) << std::endl;

	std::string line;
	getline( in, line );
	std::istringstream is( line );

	std::string tag;
	is >> tag;
	if ( !is.fail() && tag == "UPSTUB" ) {
		upstream_stub_atoms_.clear();
		while ( is >> tag ) {
			if ( tag == "DOWNSTUB" ) break;
			upstream_stub_atoms_.push_back( tag );
		}
		if ( !is.fail() && tag == "DOWNSTUB" ) {
			downstream_stub_atoms_.clear();
			while ( is >> tag ) {
				downstream_stub_atoms_.push_back( tag );
			}
		} else {
			downstream_stub_atoms_=upstream_stub_atoms_;
		}
	} else {
		set_standard_stub_atoms();
	}
	tr.Trace << "set the STUB atoms: " << std::endl;
	tr.Trace << " UPSTUB ";
	for ( AtomList::const_iterator it = upstream_stub_atoms_.begin(),
			eit = upstream_stub_atoms_.end(); it!=eit; ++it ) {
		tr.Trace << A(3, (*it));
	}
	tr.Trace << " DOWNSTUB ";
	for ( AtomList::const_iterator it = downstream_stub_atoms_.begin(),
			eit = downstream_stub_atoms_.end(); it!=eit; ++it ) {
		tr.Trace << A(3, (*it));
	}
	if ( upstream_stub_atoms_.size() != 4  || downstream_stub_atoms_.size() !=4 ) {
		utility_exit_with_message( "[ERROR] reading JumpSRFD... need 4 or NONE of the stub atoms...");
	}
	tr.Trace << std::endl;
}

} //fragment
} //core
