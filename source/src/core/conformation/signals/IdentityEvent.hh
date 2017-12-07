// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/conformation/signals/IdentityEvent.hh
/// @brief  signal a residue identity change in a Conformation (e.g. residue identity change)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_signals_IdentityEvent_hh
#define INCLUDED_core_conformation_signals_IdentityEvent_hh


// type headers
#include <core/types.hh>

// unit headers
#include <core/conformation/signals/IdentityEvent.fwd.hh>

// package headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>


namespace core {
namespace conformation {
namespace signals {


/// @brief signals a change in residue identity in a Conformation
/// @remarks When accessing residue information, take care as to which
///  data member you choose.  For almost all situations the ResidueCAP
///  'residue' should be used instead of the Conformation.  See remarks
///  below.
struct IdentityEvent : public GeneralEvent {


	// typedefs
	typedef core::Size Size;
	typedef GeneralEvent Super;


	/// @brief the type of length change
	enum Tag {
		EMPTY, // empty event, e.g. for default constructor
		INVALIDATE, // for safety, e.g. if during copy assignment something radical occurs
		RESIDUE
	};


	/// @brief default constructor
	inline
	IdentityEvent() :
		Super(),
		tag( EMPTY ),
		position( 0 ),
		residue()
	{}


	/// @brief constructor
	/// @param pos residue position
	inline
	IdentityEvent(
		Conformation const * conf,
		Tag const t,
		Size const & pos,
		Residue const * res
	) :
		Super( conf ),
		tag( t ),
		position( pos ),
		residue( res )
	{}


	/// @brief copy constructor
	inline
	IdentityEvent( IdentityEvent const & rval ) :
		Super( rval ),
		tag( rval.tag ),
		position( rval.position ),
		residue( rval.residue )
	{}


	/// @brief default destructor
	inline
	virtual
	~IdentityEvent() {}


	/// @brief copy assignment
	inline
	IdentityEvent &
	operator =( IdentityEvent const & rval ) {
		if ( this != &rval ) {
			Super::operator =( rval );
			tag = rval.tag;
			position = rval.position;
			residue = rval.residue;
		}
		return *this;
	}


	/// @brief tag indicating type of identity change
	Tag tag;

	/// @brief residue position
	Size position;

	/// @brief direct access to residue
	/// @remarks Almost always want to use this to access the residue instead of
	///  the conformation.  Calling Conformation::residue() can cause an internal
	///  update/re-sync inside Pose, which may have consequences if you're depending
	///  upon multiple residue operations to be setup (such as bond angle/length
	///  changes) prior to an internal update.
	Residue const * residue;


};


} // namespace signals
} // namespace conformation
} // namespace core


#endif /* INCLUDED_core_conformation_signals_IdentityEvent_HH */
