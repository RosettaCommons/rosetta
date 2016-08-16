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

/// @file   core/conformation/signals/GeneralEvent.hh
/// @brief  signal for a general change in a Conformation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_signals_GeneralEvent_hh
#define INCLUDED_core_conformation_signals_GeneralEvent_hh


// unit headers
#include <core/conformation/signals/GeneralEvent.fwd.hh>

// package headers
#include <core/conformation/Conformation.fwd.hh>


namespace core {
namespace conformation {
namespace signals {


/// @brief signals a general change in a Conformation
struct GeneralEvent {


	// Try to keep this class tag-less. For specific types of events
	// derive from this class.


	// typedefs
	typedef core::conformation::Conformation Conformation;


	/// @brief default constructor
	inline
	GeneralEvent() :
		conformation( NULL )
	{}


	/// @brief constructor
	inline
	GeneralEvent( Conformation const * conf ) :
		conformation( conf )
	{}


	/// @brief copy constructor
	inline
	GeneralEvent( GeneralEvent const & rval ) :
		conformation( rval.conformation )
	{}


	/// @brief default destructor
	inline
	virtual
	~GeneralEvent() {}


	/// @brief copy assignment
	inline
	GeneralEvent &
	operator =( GeneralEvent const & rval ) {
		if ( this != &rval ) {
			conformation = rval.conformation;
		}
		return *this;
	}


	/// @brief the Conformation firing the signal
	Conformation const * conformation;


};


} // namespace signals
} // namespace conformation
} // namespace core


#endif /* INCLUDED_core_conformation_signals_GeneralEvent_HH */
