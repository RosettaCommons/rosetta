// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/Frame.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_BaseCacheUnit_HH
#define INCLUDED_core_fragment_BaseCacheUnit_HH

// Unit Headers
//#include <core/fragment/BaseCacheUnit.fwd.hh>

// Package Headers

// Project Headers

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/pointer/access_ptr.fwd.hh>


// C++ STL Headers

namespace core {
namespace fragment {

class BaseCacheUnit;
typedef utility::pointer::shared_ptr< BaseCacheUnit > BaseCacheUnitOP;
typedef utility::pointer::weak_ptr< BaseCacheUnit > BaseCacheUnitAP;

class BaseCacheUnit : public utility::pointer::ReferenceCount {
public:
	virtual BaseCacheUnitOP clone() const = 0;
	virtual void remap_value( BaseCacheUnit const& source, Size source_id, Size new_id ) = 0;

	// one could call this when frags are added... that might avoid some breaking via out-of-range errors,
	// but the values would still be nonsense ... maybe out-of-range errors are preferrable
	//	virtual void register_frag_id( Size frag_id ) = 0;
};


} //fragment
} //core
#endif
