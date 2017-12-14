// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <core/types.hh>
#include <protocols/moves/ResId.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace protocols {
namespace moves {

ResId::ResId() {
	modifiable( true );
}

ResId::~ResId() = default;

ResId::ResId( core::Size const r ){
	set_resid( r );
	modifiable( true );
}

core::Size
ResId::get_resid() const{
	return( resid_ );
}

core::Size &
ResId::get_resid(){
	return( resid_ );
}

bool
ResId::modifiable() const{
	return( modifiable_ );
}

void
ResId::modifiable( bool const u ) {
	modifiable_ = u;
}

void
ResId::set_resid( core::Size const r ){
	resid_ = r;
}

/// @details a recursive function that sets the resid of obj to the resid parameter.
/// If the type of obj is a CompoundFilter or a DockDesignMover then each of the
/// members of these containers are probed to see whether they are ResId types.
/// Recursion ensures that any nesting structure would be supported.
/// Non-ResId objects will pass through this function without being changed.
void
modify_ResId_based_object( utility::pointer::ReferenceCountOP const obj, core::Size const resid ){
	using namespace protocols::filters;

	auto * resid1 = dynamic_cast< ResId * >( obj.get() );
	bool const is_this_a_resid( resid1 );
	if ( is_this_a_resid ) {
		if ( resid1->modifiable() ) {
			resid1->set_resid( resid );
		}
	}
}

} // moves
} // protocols
