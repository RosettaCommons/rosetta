// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragID.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_FragID_HH
#define INCLUDED_core_fragment_FragID_HH

// Unit Headers
#include <core/fragment/FragID.fwd.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>
#ifdef WIN32
#include <core/fragment/Frame.hh> // WIN32 INCLUDE
#endif
#include <core/fragment/FragData.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// a very lightweight class ( only two memory cells )
// copy by value, don't use OPs -- not derived from ReferenceCount.
// FragID identifies a physical fragment by its Frame and intra_frame reference number.

namespace core {
namespace fragment {

class FragID {
public:
	FragID();
	FragID( FrameCOP frame, Size frag_id );
	FragID( FragID const & );
	~FragID();
	FragID & operator = ( FragID const & rhs );

	FrameCOP frame_ptr() const;
	Frame const& frame() const;
	//Frame& frame();
	Size id() const;
	FragData const & fragment() const;
	FragDataCOP fragment_ptr() const;
	Size apply( kinematics::MoveMap const& mm, pose::Pose& pose) const;
	Size apply( pose::Pose& pose) const;

	// if we enable id != nr frame will need a map that can do frag_id --> nr
	bool is_valid() const;
	Size apply_ss( kinematics::MoveMap const& mm, std::string& ss ) const;

	/// @brief Equality operator (previously free when FragID derived from std::pair).
	bool operator == ( FragID const & other ) const;
	/// @brief Comparison operator (previously free when FragID derived from std::pair).
	bool operator <  ( FragID const & other ) const;

	/// TEMP!  soon to be private!
public:
	/// Pretend this is a std pair
	FrameCOP first;
	Size    second;

};

} //fragment
} //core

#endif
