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

/// @file   core/conformation/signals/XYZEvent.hh
/// @brief  signal for a change in XYZ coordinates in a Conformation
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_signals_XYZEvent_hh
#define INCLUDED_core_conformation_signals_XYZEvent_hh


// unit headers
#include <core/conformation/signals/XYZEvent.fwd.hh>

// package headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>


namespace core {
namespace conformation {
namespace signals {


/// @brief signals a change in XYZ coordinates in a Conformation
/// @remarks This signal should be sent when XYZs and DOFs change in
///  the atom tree inside the Conformation.  Treat the atom tree as
///  the primary subject, as the residues container is designed to
///  be update as necessary and will lag behind.
struct XYZEvent : public GeneralEvent {


	// typedefs
	typedef GeneralEvent Super;


	/// @brief default constructor
	inline
	XYZEvent() :
		Super()
	{}


	/// @brief constructor
	inline
	XYZEvent( Conformation const * conf ) :
		Super( conf )
	{}


	/// @brief copy constructor
	inline
	XYZEvent( XYZEvent const & rval ) :
		Super( rval )
	{}


	/// @brief copy assignment
	inline
	XYZEvent &
	operator =( XYZEvent const & rval ) {
		if ( this != &rval ) {
			Super::operator =( rval );
		}
		return *this;
	}


	/// @brief default destructor
	inline
	virtual
	~XYZEvent() {}


};


} // namespace signals
} // namespace conformation
} // namespace core


#endif /* INCLUDED_core_conformation_signals_XYZEvent_HH */
