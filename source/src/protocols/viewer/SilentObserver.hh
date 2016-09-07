// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_viewer_SilentObserver_hh
#define INCLUDED_protocols_viewer_SilentObserver_hh


// Unit headers
#include <protocols/viewer/SilentObserver.fwd.hh>

// Package headers
#include <core/io/silent/SilentFileData.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/Link.hh>

#include <utility/vector1.hh>
#include <iostream>


namespace protocols {
namespace viewer {

// base class
class SilentObserver : public utility::pointer::ReferenceCount {

private:

	typedef utility::pointer::ReferenceCount Super;

public:

	typedef core::pose::Pose Pose;

	/// @brief default constructor
	SilentObserver();

	/// @brief constructor
	SilentObserver( std::string const & name, bool fullatom );

	/// @brief default destructor
	~SilentObserver() override;

private: // disallow copy

	/// @brief disallow copy constructor
	// NOTE: if implementing copy constructor, remember to abstain from copying
	//       energy_event_link_ as there is no transferral of subject Pose upon copy
	SilentObserver( SilentObserver const & rval );

	/// @brief disallow copy assignment
	// NOTE: if implementing copy assignment, remember to leave energy_event_link_
	//       untouched as any current subject Pose is kept on copy assign
	SilentObserver & operator =( SilentObserver const & rval );

public :

	/// @brief attach to a Pose
	void
	attach_to( Pose & pose );

	/// @grief detach from Pose
	void
	detach_from();

	/// @brief upon receiving an EnergyEvent write to silent file
	void
	on_energy_change( core::pose::signals::EnergyEvent const & event );

private:
	int frame_count_;
	bool fullatom_;
	std::string silent_file_name_;
	core::io::silent::SilentFileDataOP sfd_;
	utility::signals::Link energy_event_link_;

};

} // viewer
} // protocols


#endif
