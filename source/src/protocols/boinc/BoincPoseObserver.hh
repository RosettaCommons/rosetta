// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/boinc/BoincPoseObserver.hh
/// @brief  BoincPoseObserver class
/// @author Will Sheffler, David E Kim


#ifndef INCLUDED_protocols_boinc_BoincPoseObserver_hh
#define INCLUDED_protocols_boinc_BoincPoseObserver_hh


// Unit headers
#include <protocols/boinc/BoincPoseObserver.fwd.hh>

#include <protocols/boinc/boinc_shmem.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>

#include <utility/signals/Link.hh>

// C++ Headers


namespace protocols {
namespace boinc {

// base class
class BoincCurrentPoseObserver : public utility::pointer::ReferenceCount {

public: // typedefs

	typedef core::pose::Pose Pose;
	typedef core::pose::signals::ConformationEvent ConformationEvent;

public:

	/// @brief default constructor
	BoincCurrentPoseObserver();

	/// @brief default destructor
	virtual ~BoincCurrentPoseObserver();

private: // disallow copy

	/// @brief disallow copy constructor
	// NOTE: if implementing copy constructor, do not copy the 'conf_event_link_'
	//       as there is no transferal of subject Pose on copy construct
	BoincCurrentPoseObserver( BoincCurrentPoseObserver const & rval );

	/// @brief disallow copy assignment
	// NOTE: if ConformationViewer copy assignment, remember to leave 'conf_event_link_'
	//       untouched as any current subject Pose is kept on copy assign
	BoincCurrentPoseObserver & operator =( BoincCurrentPoseObserver const & rval );

public:

	/// @brief attach to Pose
	void
	attach_to( Pose & pose );

	/// @brief detach from Pose
	void
	detach_from();

	/// @brief on receiving ConformationEvent, copy Pose to boinc shared memory
	void
	on_conf_change( ConformationEvent const & event );

private:

	BoincSharedMemory * shmem_;

	utility::signals::Link conf_event_link_;

};

} // namespace boinc
} // namespace protocols


#endif // INCLUDED_protocols_boinc_BoincPoseObserver_HH
