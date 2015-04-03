// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/ResId.hh
/// @brief A simple base class providing common access to a residue identity for movers and filters.
///        The anticipated usage is to derive movers and filters multiply, from the relevant base class and
///        from ResId and make these movers/filters use the residue id encapsulated in this class.
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_moves_ResId_hh
#define INCLUDED_protocols_moves_ResId_hh

#include <core/types.hh>
#include <protocols/moves/ResId.fwd.hh>

#include <utility/pointer/ReferenceCount.fwd.hh>


namespace protocols {
namespace moves {


//// Override the get and set methods in containers of ResId objects (see CompoundStatementFilter and DockDesignMover for examples)
class ResId{
public:
	ResId();
	ResId( core::Size const r );
	virtual void set_resid( core::Size const r );
	virtual core::Size get_resid() const;
	virtual core::Size & get_resid();
	/// @brief should another method be able to modify resid_. This is used by modify_ResId_based_object as a test
	virtual bool modifiable() const;
	virtual void modifiable( bool const u );
	virtual ~ResId();
private:
	core::Size resid_;
	bool modifiable_;// defaults to true
};

/// @brief Checks whether a referencecount object is a derived from ResId and if so, sets its resid
void
modify_ResId_based_object( utility::pointer::ReferenceCountOP const obj, core::Size const resid );

} // moves
} // protocols

#endif
