// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/PoseMetricContainer.hh
/// @brief  container for managing PoseMetricCalculators
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_PoseMetricContainer_hh
#define INCLUDED_core_pose_metrics_PoseMetricContainer_hh


// unit headers
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>

// package headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <basic/MetricValue.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
#include <map>
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace metrics {


/// @brief container for managing PoseMetricCalculators
/// @details Attaches to a Pose and listens to Pose signals.
///  It then notifies stored PoseMetricCalculators to recalculate as-needed.
class PoseMetricContainer : public utility::pointer::ReferenceCount {


private: // typedefs

	typedef std::map< std::string, PoseMetricCalculatorOP > Name2Calculator;


public: // construct/destruct


	/// @brief default destructor
	PoseMetricContainer();


	/// @brief copy constructor
	/// @warning observation of subject Pose in src not duplicated
	PoseMetricContainer( PoseMetricContainer const & src );


	/// @brief default destructor
	virtual ~PoseMetricContainer();


	/// @brief copy assignment
	/// @warning container keeps observing whatever Pose it was originally observing
	PoseMetricContainer & operator =( PoseMetricContainer const & src );


public: // calculator management


	/// @brief clear the list of metric calculators, reset flags
	void clear();


	/// @brief get a value out of a PoseMetricCalculator
	// TODO: should we really be passing a Pose pointer in...?  Consider using the internal pose_ptr_;
	void get( std::string const & calculator_name, std::string const & key, basic::MetricValueBase & val, Pose const & this_pose );


	/// @brief return a string with the results of a PoseMetricCalculator
	std::string print( std::string const & calculator_name, std::string const & key, Pose const & this_pose );


public: // observer interface


	/// @brief attach to a Pose
	void attach_to( core::pose::Pose & pose );


	/// @brief detach from Pose
	void detach_from();


	/// @brief is observing a Pose?
	/// @return the Pose if observing, otherwise NULL
	Pose const * is_observing();


	/// @brief upon receiving a pose::signals::DestructionEvent, detaches
	void on_pose_destruction_change( core::pose::signals::DestructionEvent const & event );


	/// @brief upon receiving a pose:signals::EnergyEvent, sets flag telling
	///  calculators to refresh energy based calculations
	void on_pose_energy_change( core::pose::signals::EnergyEvent const & event );


	/// @brief upon receiving pose::signals::ConformationEvent, sets flag telling
	///  calculators to refresh structure based calculations
	void on_pose_conf_change( core::pose::signals::ConformationEvent const & event );


private: // calculator update


	/// @brief set PoseMetricCalculator structure changed flags
	void process_structure_change();


	/// @brief set PoseMetricCalculator energy changed flags
	void process_energy_change();


	/// @brief get a PoseMetricCalculator by name
	PoseMetricCalculatorOP get_calculator( std::string const & calculator_name );


	/// @brief add a PoseMetricCalculator by name
	void add_calculator( std::string const & calculator_name );


	/// @brief clone calculators, primarily for copy
	void clone_calculators( Name2Calculator const & calcs );


private: // data


	/// @brief pointer to the Pose being watched
	core::pose::Pose const * pose_ptr_;


	/// @brief flag for structure change
	bool structure_is_outdated_;


	/// @brief flag for energy change
	bool energies_are_outdated_;


	/// @brief the list of metric calculators
	Name2Calculator metric_calculators_;

};


} // namespace metrics
} // namespace pose
} // namespace core


#endif
