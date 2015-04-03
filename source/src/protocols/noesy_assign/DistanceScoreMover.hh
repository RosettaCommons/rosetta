// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_DistanceScoreMover_HH
#define INCLUDED_protocols_noesy_assign_DistanceScoreMover_HH


// Unit Headers
#include <protocols/noesy_assign/DistanceScoreMover.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
//#include <protocols/noesy_assign/PeakAssignment.hh>
#include <core/scoring/constraints/Constraint.hh>
// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/noesy_assign/CrossPeakList.fwd.hh>

#ifdef WIN32
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#endif

namespace protocols {
namespace noesy_assign {

/// @brief maintains a list of constraints_ (each PeakAssignment yields one) and peak_constraints_ ( each
/// cross peak with multiple assignments creates one AmbiguousNMRConstraint ).
/// the DistanceMover (prepare scoring and apply directly interacts with the Dk term in CrossPeak (set_decoy_compatibility)


class DistanceScoreMover : public protocols::moves::Mover {
public:
	typedef utility::vector1< core::pose::PoseOP > PoseVector;

  DistanceScoreMover( CrossPeakList&, core::pose::Pose const& pose, core::Real dcut );

	/// @brief set decoy_compatibility in PeakAssignments to zero
	///
  void prepare_scoring( bool use_for_calibration = false );

	/// @brief sum up decoy_compatibility score in PeakAssignments
  void apply( core::pose::Pose& pose );

	//	void find_violators_with_individual_dist_cutoff( PoseVector poses );
	/// @brief normalize decoy_compatibility of PeakAssignments by count_decoys_
  void finalize_scoring() const;

	//core::Real compute_violation_percentage() const;
	//void eliminate_violated_constraints() const;

	virtual std::string get_name() const { return "DistanceScoreMover"; }

	void set_dcut( core::Real setting ) {
		dcut_ = setting;
	}

private:
  CrossPeakList& cross_peaks_;
  core::Size count_decoys_; //how many decoys for scoring
  core::Size nr_assignments_;
  typedef utility::vector1< core::scoring::constraints::ConstraintOP > SingleConstraints;
  SingleConstraints constraints_;


	//use this if we have a distance cutoff for each individual peak, based on the structural variation at that point.
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool individual_distance_cutoff_;
	typedef utility::vector1< core::Real > DistanceBoundVector;
	typedef utility::vector1< DistanceBoundVector > AllStructureDistanceBoundVector;
	AllStructureDistanceBoundVector all_dist_buf_;
	//	typedef utility::vector1< core::scoring::constraints::AmbiguousNMRConstraintOP > PeakConstraints;
	//	PeakConstraints peak_constraints_;


	/// @brief count for each peak how many decoys have violated
	typedef utility::vector1< core::Size > PeakViolationCounts;
	PeakViolationCounts peak_violation_counts_;

	/// @brief cumulative sum of peak_violation_counts_;
	core::Size total_violation_count_;

	//	typedef utility::vector1< core::Real > VectorReal;
	core::Real total_assigned_distances_;


	core::Size active_peaks_; //count nr peaks that have assignment with sufficient Vk

  core::Real final_dist_power_; //eta in Cyana-paper

	bool used_for_calibration_; //if used_for_calibration --> filter peaks and do not update decoy_compatibility score
	core::Real dcut_;
};

}
}

#endif

