// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/conformation/signals/LengthEvent.hh
/// @brief  signal for a change in length of residues in a Conformation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Rocco Moretti (rmorettiase@gmail.com)

// unit headers
#include <core/conformation/signals/LengthEvent.hh>

// package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/Conformation.hh>

// utility headers
#include <utility/exit.hh>

#include <iostream>

namespace core {
namespace conformation {
namespace signals {

/// @brief constructor
/// @param t type of length change
/// @param pos residue position
LengthEvent::LengthEvent(
	Conformation const * conf,
	Tag const t,
	Size const & pos,
	int const & len_chg,
	Residue const * res
) :
	Super( conf ),
	tag( t ),
	position( pos ),
	conformation_size( 0 ),
	length_change( len_chg ),
	residue( res )
{
	if ( conf ) {
		conformation_size = conf->size();
	}
#ifndef NDEBUG
	check_consistency();
#endif
}

void
LengthEvent::check_consistency() const {
	switch( tag ) {
	case RESIDUE_APPEND:
	case RESIDUE_PREPEND :
		runtime_assert( length_change >= 0 );
		break;
	case RESIDUE_DELETE :
		runtime_assert( length_change <= 0 );
		break;
	default :
		break; // In general, there's no need of this.
	}
}

std::ostream &
operator<<( std::ostream & out, LengthEvent const & event ) {
	out << "LengthEvent of type ";
	switch( event.tag ) {
	case LengthEvent::EMPTY :
		out << "EMPTY";
		break;
	case LengthEvent::INVALIDATE :
		out << "INVALIDATE";
		break;
	case LengthEvent::RESIDUE_APPEND :
		out << "APPEND";
		break;
	case LengthEvent::RESIDUE_PREPEND :
		out << "PREPEND";
		break;
	case LengthEvent::RESIDUE_DELETE :
		out << "DELETE";
		break;
	default :
		out << "UNKNOWN";
		break;
	}
	out << " at " << event.position << " of length " << event.length_change;
	out << " on conformation " << event.conformation << " (final size " << event.conformation_size << " ) and residue " << event.residue;
	return out;
}

} // namespace signals
} // namespace conformation
} // namespace core

