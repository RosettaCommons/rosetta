// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/Frame.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/JumpingFrame.hh>

// Package Headers

// Project Headers
#include <core/id/SequenceMapping.hh>
//#include <core/chemical/ChemicalManager.hh>
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/conformation/ResidueFactory.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


// Utility headers
#include <utility/excn/Exceptions.hh>


namespace core {
namespace fragment {

using namespace ObjexxFCL::format;

// FragDataOP generate_fragdata( SingleResidueFragDataOP frag_res_type, SingleResidueFragDataOP jump_frag_type  ) {
//   if ( nr_frags() ) return frag_data_.front()->clone();
//   Size const length( bWithTorsion ? 3 : 1 );
//   if ( bWithTorsion ) {
//     BBTorsionSRFDOP start =  new BBTorsionSRFD( 3, 'E', 'X' );
//     start->set_torsion( 1, it->phi( iStart ) );
//     start->set_torsion( 2, it->psi( iStart ) );
//     start->set_torsion( 3, it->omega( iStart ) );

//     frags.back()->add_residue( start );
//   }
//   frags.back()->add_residue( new JumpSRFD( it->rt_, it->atoms_downstream_, it->atoms_upstream_, 'X' ) );
//   if ( bWithTorsion ) {
//     BBTorsionSRFDOP stop =  new BBTorsionSRFD( 3, 'E', 'X' );
//     stop->set_torsion( 1, it->phi( iStop ) );
//     stop->set_torsion( 2, it->psi( iStop ) );
//     stop->set_torsion( 3, it->omega( iStop ) );

//     frags.back()->add_residue( stop );
//   }

// }
bool
NonContinuousFrame::align( core::id::SequenceMapping const& map ) {
	bool success = Parent::align( map );
	for ( auto it = pos_.begin(),
			eit = pos_.end(); it!=eit && success; ++it )  {
		Size newpos( map[ *it ] );
		if ( newpos > 0 ) {
			*it = newpos;
		} else return false;
	}
	return success;
}

void NonContinuousFrame::show( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	out << type() << " ";
	show_pos( out );
	out << std::endl;
	show_fragments( out );
}

void NonContinuousFrame::read( std::istream &in ) {
	using namespace ObjexxFCL::format;
	Size pos;
	pos_.clear();
	std::string line;
	getline( in, line );
	std::istringstream line_stream( line );
	while ( line_stream >> pos )  {
		pos_.push_back( pos );
	}
	init_length( pos_.front(), pos_.back(), pos_.size() );
}

void NonContinuousFrame::show_pos( std::ostream &out ) const {
	for ( core::Size po : pos_ ) {
		out << RJ( 3, po ) << " ";
	}
}

void NonContinuousFrame::shift_to( core::Size setting ) {
	using core::Size;
	int offset = setting-pos_.front();
	shift_by( offset );
}

void NonContinuousFrame::shift_by( int offset ) {
	Parent::shift_by( offset );
	for ( core::Size & po : pos_ ) {
		int new_pos = po + offset;
		if ( new_pos < 1 ) {
			std::ostringstream msg;
			msg << "offset " << offset << " would shift at least one position of Frame " << *this << " to negative or zero" << std::endl;
			throw utility::excn::EXCN_RangeError( msg.str() );
		}
		po = new_pos;
	}
}


//  --- cloning with frags is taken care of correctly by base-class
// it is enough that clone() itself is virtually overloaded
/// @brief clone method, new frame with same alignment position, fragments are not copied!
FrameOP NonContinuousFrame::clone() const {
	return FrameOP( new NonContinuousFrame( start(), end(), length() ) );
}

// /// @brief clone method, new frame with same alignment position, fragments are not copied!
// NonContinuousFrameOP NonContinuousFrame::clone_with_frags() const {
//  NonContinuousFrameOP newFrame = clone();// new Frame( start(), end(), length() );
//  *newFrame = *this; //usually that is enough
//  return newFrame;
// }

// /// @brief clone method, new frame with same alignment position, one fragments is copied as template ( valid() == false )
// NonContinuousFrameOP NonContinuousFrame::clone_with_template() {
//  NonContinuousFrameOP newFrame = clone();// new Frame( start(), end(), length() );
//  if ( nr_frags() ) {
//   newFrame->frag_list_.push_back( frag_list_[ 1 ]->clone() );
//   newFrame->frag_list_[ 1 ]->set_valid( false );
//  }
//  return newFrame;
// }


}
}
