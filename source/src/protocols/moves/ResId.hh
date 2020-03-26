// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/ResId.hh
/// @brief A simple base class providing common access to a residue identity for movers and filters.
///        The anticipated usage is to derive movers and filters multiply, from the relevant base class and
///        from ResId and make these movers/filters use the residue id encapsulated in this class.
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_moves_ResId_hh
#define INCLUDED_protocols_moves_ResId_hh

#include <core/types.hh>
#include <protocols/moves/ResId.fwd.hh>
#include <core/pose/ResidueIndexDescription.fwd.hh>
#include <utility/VirtualBase.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace moves {


//// Override the get and set methods in containers of ResId objects (see CompoundStatementFilter and DockDesignMover for examples)
class ResId{
public:
	ResId() = default;
	ResId( core::Size const r );
	virtual void set_resid( core::Size const r );
	virtual void set_resid( core::pose::ResidueIndexDescriptionCOP r );
	virtual core::Size get_resid( core::pose::Pose const & ) const;
	/// @brief should another method be able to modify resid_. This is used by modify_ResId_based_object as a test
	bool modifiable() const;
	void modifiable( bool const u );
private:
	core::pose::ResidueIndexDescriptionCOP resid_;
	bool modifiable_ = true;
};

/// @brief Checks whether a referencecount object is a derived from ResId and if so, sets its resid
void
modify_ResId_based_object( utility::VirtualBaseOP const obj, core::Size const resid );

/// @brief Checks whether a referencecount object is a derived from ResId and if so, sets its resid
void
modify_ResId_based_object( utility::VirtualBaseOP const obj, core::pose::ResidueIndexDescriptionCOP r );

} // moves
} // protocols

#endif
