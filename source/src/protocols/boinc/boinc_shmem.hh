// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/boinc/boinci_shmem.hh
/// @brief  Boinc shared memory structure for graphics
/// @author David Kim (dekim@u.washington.edu)


#ifndef INCLUDED_protocols_boinc_boinc_shmem_hh
#define INCLUDED_protocols_boinc_boinc_shmem_hh


#include <core/types.hh>

#ifdef BOINC
#include "boinc_api.h"
#endif

namespace protocols {
namespace boinc {

const std::size_t MAX_NATIVE_POSE_RESIDUES = 500; // large natives cause sluggish graphics
const std::size_t MAX_SYMM_POSE_RESIDUES = 500;  // only the asymmetric unit will be displayed if greater than this
const std::size_t POSE_BUFSIZE = 5000000; // size should depend on the largest pose expected to be used //1000000; //99999;
const std::size_t TEXT_BUFSIZE = 255; //99999;
const std::size_t WU_DESC_TEXT_BUFSIZE = 1024; //99999;

struct BoincSharedMemory {
	double update_time;
	double fraction_done;
	double cpu_time;

#ifdef BOINC
	BOINC_STATUS status;
#endif

// Lets save the info necessary to reproduce rosetta++ graphics
// rmsds will be calculated by the graphics app against the native_pose below

	// current pose
	int current_pose_exists;
	char current_pose_buf[POSE_BUFSIZE];
	
	// current pose ghost
	int current_pose_ghost_exists;
	char current_pose_ghost_buf[POSE_BUFSIZE];

	// accepted
	int last_accepted_pose_exists;
	char last_accepted_pose_buf[POSE_BUFSIZE];

	// low energy
	int low_energy_pose_exists;
	char low_energy_pose_buf[POSE_BUFSIZE];

	// native pose
	char native_pose_buf[POSE_BUFSIZE];
	int native_pose_exists;

	// monte carlo total step count
	unsigned int total_mc_trial_count;

  // scores
	core::Real low_energy;
	core::Real last_accepted_energy;

	unsigned int low_energy_update_cnt;

	// model nstruct
	unsigned int model_count;
	// model low energy
	float model_low_energy;
	// model rmsd
	float model_low_energy_rmsd;

	// job info
	char job_type_text[TEXT_BUFSIZE];

	// stage info
	char mover_type_text[TEXT_BUFSIZE];

	// work unit description
	char wu_desc_buf[WU_DESC_TEXT_BUFSIZE];
	int wu_desc_exists;

	// monte carlo mover step count
	unsigned int mover_mc_trial_count;

	// Should we randomly cycle appearance?
	bool randomly_cycle_appearance;
	
	// Has the main app initialized the shared memory object completely, so that the graphics app can read what it needs?
	bool fully_initialized;

};


} // namespace boinc
} // namespace protocols

#endif

