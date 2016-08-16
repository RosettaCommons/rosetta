// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/boinc/BoincPoseObserver.cc
/// @brief  BoincPoseObserver class

#ifdef BOINC_GRAPHICS

// Unit headers
#include <protocols/boinc/BoincPoseObserver.hh>

#include <protocols/boinc/boinc.hh>

#include <core/pose/Pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/serialization/serialize_pose.hh>


namespace protocols {
namespace boinc {

/// @brief default constructor
BoincCurrentPoseObserver::BoincCurrentPoseObserver()
{
	shmem_ = protocols::boinc::Boinc::get_shmem();
	is_ghost_ = false;
}

/// @brief Constructor for observers for drawing a "ghost" overlaid on the current pose.
BoincCurrentPoseObserver::BoincCurrentPoseObserver( bool const is_ghost )
{
	shmem_ = protocols::boinc::Boinc::get_shmem();
	is_ghost_ = is_ghost;
}

/// @brief default destructor
BoincCurrentPoseObserver::~BoincCurrentPoseObserver()
{
	detach_from();
}


/// @brief attach to Pose
void
BoincCurrentPoseObserver::attach_to( Pose & pose )
{
	detach_from();
	conf_event_link_ = pose.attach_conformation_obs( &BoincCurrentPoseObserver::on_conf_change, this );
}


/// @brief detach from Pose
void
BoincCurrentPoseObserver::detach_from()
{
	conf_event_link_.invalidate();
}


/// @brief on receiving an ConformationEvent, copy Pose to boinc shared memory
void
BoincCurrentPoseObserver::on_conf_change(
	ConformationEvent const & event
)
{
	static int count = 0;
	if ( count > 0 ) {
		count--;
		return;
	}
	count = SKIP_FOR_EFFICIENCY;
	if ( !shmem_ ) return;
	if ( !Boinc::trywait_semaphore() ) {
		boinc_begin_critical_section();
		if ( event.pose->total_residue() > 0 ) {
			core::io::serialization::BUFFER b((char*)( is_ghost_ ? &shmem_->current_pose_ghost_buf : &shmem_->current_pose_buf ),POSE_BUFSIZE);
			if ( core::pose::symmetry::is_symmetric( *event.pose ) && event.pose->total_residue() > MAX_SYMM_POSE_RESIDUES ) {
				core::pose::Pose pose;
				core::pose::symmetry::extract_asymmetric_unit(*event.pose, pose);
				write_binary(pose,b);
			} else {
				write_binary(*event.pose,b);
			}
			if ( !is_ghost_ ) {
				shmem_->current_pose_exists = 1;
			} else {
				shmem_->current_pose_ghost_exists = 1;
			}
		}
		boinc_end_critical_section();
		Boinc::unlock_semaphore();
	}
}


} // namespace boinc
} // namespace protocols

#endif
