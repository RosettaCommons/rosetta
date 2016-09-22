// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PeakAssignmentParametersList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakAssignmentParameters_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignmentParameters_HH


// Unit Header

// Package Headers

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

//// C++ headers
#include <iostream>
#include <iosfwd>
#include <string>

#ifdef MULTI_THREADED
#include <atomic>
#include <mutex>
#endif

namespace protocols {
namespace noesy_assign {

// THIS SHOULD NOT BE A SINGLETON! ACCESS TO THIS CLASS IS *HIGHLY* NON-THREADSAFE
class PeakAssignmentParameters { //: public utility::pointer::ReferenceCount {

private:
	static bool options_registered_;
	PeakAssignmentParameters() {}; //private constructor

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PeakAssignmentParameters * create_singleton_instance();

public:
	static void register_options();
	static void set_cycle( core::Size );
	void show( std::ostream& ) const;
	void show_on_tracer() const;

	// KILL THESE - KILL THEM WITH FIRE!
	static PeakAssignmentParameters const* get_instance();
	static PeakAssignmentParameters * get_nonconst_instance();

	static void reset();
private:
	void set_options_from_cmdline( core::Size cycle = 0 );
	/// Singleton instance pointer
#if defined MULTI_THREADED
	static std::atomic< PeakAssignmentParameters * > instance_;
#else
	static PeakAssignmentParameters * instance_;
#endif

	core::Size cycle_selector_;


#ifdef MULTI_THREADED
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif

public:
	/* maybe make all options const
	make cmd-line options vectors for 7 cycles and have static function that
	switches cycles by throwing away one instance and creating a new one with the given
	cycle number and setting all const-params in the constructor from the cmd-line.

	all good, unless I want to change the things from somewhere else than cmdline...
	*/

	//cycle independent
	core::Real chemshift_overlap_weight_; //Gamma, eq. (4)
	bool ignore_resonancefile_tolerances_;
	bool ignore_resonancefile_intensities_;
	//  core::Real dmax_; //unused
	core::Real vmin_;
	core::Real vmax_;
	core::Real nmax_;
	core::Real nr_conformers_violatable_; //Mvio in fraction of nr_conformers
	core::Real network_reswise_high_;
	core::Real centroid_mapping_distance_padding_;
	//cycle dependent
	core::Real symmetry_compliance_weight_; //T, eq. (5)
	core::Real covalent_compliance_weight_; //O
	//obsolet  core::Real decoy_compatibility_exponent_; //eta, eq. (6)

	core::Real min_contribution_symmetric_peaks_;
	core::Real smax_;
	core::Real dcut_;
	core::Real dcalibrate_;

	bool use_local_distviol_;
	core::Real local_distviol_range_; //how many (in percent) decoys at both ends of the range are ignored to calculate max_extension
	core::Real local_distviol_global_buffer_;
	core::Real local_distviol_global_factor_;
	core::Real local_distviol_cutoff_;
	core::Real local_distviol_cutoff_buffer_;

	core::Real network_reswise_min_; //N_bar_min
	core::Real network_atom_min_; //N_min per atom
	core::Real calibration_target_;
	bool calibration_ignore_eliminated_peaks_;
	bool atom_dependent_calibration_;
	core::Real min_volume_; //minimum volume contribution

	core::Real cst_strength_; //for  1/cst_strength ->sigma for BoundFunc

	bool no_network_;
	bool network_use_all_covalent_atoms_;
	bool network_include_reverse_dir_;
	bool network_allow_same_residue_connect_;
	std::string network_mode_;
	bool map_to_cen_atom_;

	core::Real calibration_convergence_;
	core::Real calibration_max_noe_dist_;
	core::Real calibration_stop_nudging_;
	core::Real calibration_start_nudging_;
	core::Real calibration_max_nudging_;
	bool calibration_eliminate_;
	bool calibration_use_median_;
	core::Size calibration_cycles_;

	utility::vector1<core::Real> prob_sigmoid_tau_;
	utility::vector1<core::Real> prob_sigmoid_m_;
	utility::vector1<core::Real> prob_sigmoid_w_;
	utility::vector1<core::Real> prob_level_;

};

}
}

#endif
