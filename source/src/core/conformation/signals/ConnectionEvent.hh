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

/// @file   core/conformation/signals/ConnectionEvent.hh
/// @brief  signal a change in the connection with a Conformation object,
///         e.g. destruction or transfer
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_signals_ConnectionEvent_hh
#define INCLUDED_core_conformation_signals_ConnectionEvent_hh


// unit headers
#include <core/conformation/signals/ConnectionEvent.fwd.hh>

// package headers
#include <core/conformation/Conformation.fwd.hh>


namespace core {
namespace conformation {
namespace signals {


/// @brief signal a change in the connection with a Conformation object, e.g.
///  destruction or transfer
/// @remarks SUGGESTION: Try to use the Link management provided by the
///  SignalHub instead of listening to this event, as it typically makes
///  managing connections much easier.
struct ConnectionEvent { // do not derive from GeneralEvent


	// typedefs
	typedef core::conformation::Conformation Conformation;


	/// @brief the type of Conformation lifetime event
	/// @details Current tags are as follows:
	///  <ul>
	///    <li> 'EMPTY' - null event for default ConnectionEvent constructor (no point
	///         in watching for this non-event)
	///    <li> 'DISCONNECT' - force disconnect (e.g. the Conformation is getting destroyed)
	///    <li> 'TRANSFER' - the connection is getting transferred to a
	///         new Conformation, so if any observers are storing Conformation
	///         pointers (e.g. those not using Link management) then they
	///         should discard the existing Conformation pointer and swap it
	///         with the one provided by the Event.
	///  </ul>
	enum Tag {
		EMPTY, // empty event for default constructor
		DISCONNECT, // force disconnect (e.g. the Conformation is getting destroyed)
		TRANSFER // if observer stores a Conformation pointer
	};


	/// @brief default constructor
	inline
	ConnectionEvent() :
		conformation( NULL ),
		tag( EMPTY )
	{}


	/// @brief constructor
	/// @param[in] conf The Conformation firing the signal.
	/// @param[in] t The tag specifying the type of signal.
	inline
	ConnectionEvent(
		Conformation const * conf,
		Tag const t
	) :
		conformation( conf ),
		tag( t )
	{}


	/// @brief copy constructor
	inline
	ConnectionEvent( ConnectionEvent const & rval ) :
		conformation( rval.conformation ),
		tag( rval.tag )
	{}


	/// @brief default destructor
	inline
	virtual
	~ConnectionEvent() {}


	/// @brief copy assignment
	inline
	ConnectionEvent &
	operator =( ConnectionEvent const & rval ) {
		if ( this != &rval ) {
			conformation = rval.conformation;
			tag = rval.tag;
		}
		return *this;
	}


	/// @brief the Conformation firing the signal
	Conformation const * conformation;


	/// @brief tag indicating type of connection change
	Tag tag;


};


} // namespace signals
} // namespace conformation
} // namespace core


#endif /* INCLUDED_core_conformation_signals_ConnectionEvent_HH */
