// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_PoseMetricCalculatorBase_hh
#define INCLUDED_core_pose_metrics_PoseMetricCalculatorBase_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace pose {
namespace metrics {

// Note: PoseMetricCalculator is an abstract base class for
//           "StuctureDependentCalculator" and "EnergyDependentCalculator".
// Note:  "SequenceDependentCalculator" may come later...

// Other assorted calculator types derive from
// StuctureDependentCalculator and SequenceDependentCalculator
// These derived classes should hold cached data, and define the following three functions:
//   "recompute" (update all cached data)
//   "lookup" (return the value of a piece of cached data inside a MetricValue)
//   "print" (return the value of a piece of cached data as a string)
// Note: Derived classes should NOT redefine "get"

class PoseMetricCalculator : public utility::pointer::ReferenceCount {
public:

	PoseMetricCalculator() {};
	//these functions are not defined and thus should not be in the header
	//PoseMetricCalculator( PoseMetricCalculator const & src );
	//PoseMetricCalculator const & operator=( PoseMetricCalculator const & src );
	virtual PoseMetricCalculatorOP clone() const = 0;

	virtual void notify_structure_change() { return; };
	// virtual void notify_sequence_change() { return; };
	virtual void notify_energy_change() { return; };

	virtual void get( std::string const & key, basic::MetricValueBase & val, Pose const & this_pose ) = 0;

	virtual std::string get( std::string const & key, Pose const & this_pose ) = 0;

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const = 0;

	virtual std::string print( std::string const & key ) const = 0;

	virtual void recompute( Pose const & this_pose ) = 0;

};


class StructureDependentCalculator : public PoseMetricCalculator {
public:
	StructureDependentCalculator() : PoseMetricCalculator(), structure_is_outdated_(true) {};
	void notify_structure_change() { structure_is_outdated_ = true; };
	void get( std::string const & key, basic::MetricValueBase & val, Pose const & this_pose ) {
		if ( structure_is_outdated_ ) recompute( this_pose );
		structure_is_outdated_ = false;
		lookup( key, &val );
	};
	std::string get( std::string const & key, Pose const & this_pose ) {
		if ( structure_is_outdated_ ) recompute( this_pose );
		structure_is_outdated_ = false;
		return print( key );
	};
protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const = 0;
	virtual std::string print( std::string const & key ) const = 0;
	virtual void recompute( Pose const & this_pose ) = 0;
private:
	bool structure_is_outdated_;
};

/*
class SequenceDependentCalculator : public PoseMetricCalculator {
public:
SequenceDependentCalculator() : PoseMetricCalculator(), sequence_is_outdated_(true) {};
void notify_sequence_change() { sequence_is_outdated_ = true; };
void get( std::string const & key, basic::MetricValueBase & val, Pose const & this_pose ) {
if ( sequence_is_outdated_ ) recompute( this_pose );
sequence_is_outdated_ = false;
lookup( key, &val );
};
std::string get( std::string const & key, Pose const & this_pose ) {
if ( sequence_is_outdated_ ) recompute( this_pose );
sequence_is_outdated_ = false;
return print( key );
};
protected:
virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const = 0;
virtual std::string print( std::string const & key ) const = 0;
virtual void recompute( Pose const & this_pose ) = 0;
private:
bool sequence_is_outdated_;
};
*/

class EnergyDependentCalculator : public PoseMetricCalculator {
public:
	EnergyDependentCalculator() : PoseMetricCalculator(), energies_are_outdated_(true) {};
	void notify_energy_change() { energies_are_outdated_ = true; };
	void get( std::string const & key, basic::MetricValueBase & val, Pose const & this_pose ) {
		if ( energies_are_outdated_ ) recompute( this_pose );
		energies_are_outdated_ = false;
		lookup( key, &val );
	};
	std::string get( std::string const & key, Pose const & this_pose ) {
		if ( energies_are_outdated_ ) recompute( this_pose );
		energies_are_outdated_ = false;
		return print( key );
	};
protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const = 0;
	virtual std::string print( std::string const & key ) const = 0;
	virtual void recompute( Pose const & this_pose ) = 0;
private:
	bool energies_are_outdated_;
};


} // namespace metrics
} // namespace pose
} // namespace core


#endif
