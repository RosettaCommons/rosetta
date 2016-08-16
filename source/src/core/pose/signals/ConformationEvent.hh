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

/// @file   core/pose/signals/ConformationEvent.hh
/// @brief  signal for a conformation change in a Pose
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_signals_ConformationEvent_hh
#define INCLUDED_core_pose_signals_ConformationEvent_hh


// unit headers
#include <core/pose/signals/ConformationEvent.fwd.hh>

// package headers
#include <core/pose/signals/GeneralEvent.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace signals {


/// @brief signals a conformation change in a Pose
struct ConformationEvent : public GeneralEvent {


	// typedefs
	typedef GeneralEvent Super;


	/// @brief default constructor
	inline
	ConformationEvent() :
		Super()
	{}


	/// @brief constructor
	inline
	ConformationEvent( Pose const * pose ) :
		Super( pose )
	{}


	/// @brief copy constructor
	inline
	ConformationEvent( ConformationEvent const & rval ) :
		Super( rval )
	{}


	/// @brief default destructor
	inline
	virtual
	~ConformationEvent() {}


	/// @brief copy assignment
	inline
	ConformationEvent &
	operator =( ConformationEvent const & rval ) {
		if ( this != &rval ) {
			Super::operator =( rval );
		}
		return *this;
	}


};


} // namespace signals
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_signals_ConformationEvent_HH */
