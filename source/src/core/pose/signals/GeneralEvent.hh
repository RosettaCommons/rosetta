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

/// @file   core/pose/signals/GeneralEvent.hh
/// @brief  signal for a general change in a Pose
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pose_signals_GeneralEvent_hh
#define INCLUDED_core_pose_signals_GeneralEvent_hh


// unit headers
#include <core/pose/signals/GeneralEvent.fwd.hh>


// Package headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace signals {


/// @brief signals a general change in a Pose
struct GeneralEvent {


	// Try to keep this tag-less. For specific types of events
	// derive from this class.


	// typedefs
	typedef core::pose::Pose Pose;


	/// @brief default constructor
	inline
	GeneralEvent() :
		pose( NULL )
	{}


	/// @brief constructor
	inline
	GeneralEvent( Pose const * pose ) :
		pose( pose )
	{}


	/// @brief copy constructor
	inline
	GeneralEvent( GeneralEvent const & rval ) :
		pose( rval.pose )
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
			pose = rval.pose;
		}
		return *this;
	}


	/// @brief the Pose firing the signal
	Pose const * pose;

	/// @brief PyRosetta work around
	Pose const & getPose() { return *pose; };
};


} // namespace signals
} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_signals_GeneralEvent_HH */
