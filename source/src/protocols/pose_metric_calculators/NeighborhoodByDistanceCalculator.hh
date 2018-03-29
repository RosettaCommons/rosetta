// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NeighborhoodByDistanceCalculator.hh
/// @brief NeighborhoodByDistanceCalculator can determine all the neighbors of group of residues within a certain distance.
/// @author Steven Lewis

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborhoodByDistanceCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborhoodByDistanceCalculator_hh

//Unit headers
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/pose_metric_calculators/NeighborhoodByDistanceCalculator.fwd.hh>


#include <core/pose/Pose.fwd.hh>
#include <basic/MetricValue.fwd.hh>

//Utility headers
#include <core/types.hh>


//C++ headers
#include <set>

#include <utility/vector1.hh>
#include <map>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {

/// @details this calculator determines the number and resids of residues within X angstroms of a group of given residues.  Its intended purpose is the backend for a TaskOperation that allows one to construct a PackerTask based on neighborhoods around a set of particular residues.  It can return its set of central residues, the total count of their neighbors as determined by the sub-calculators (inclusive of the central residues), and the identities of those neighbors.
class NeighborhoodByDistanceCalculator : public core::pose::metrics::StructureDependentCalculator {

public:
	typedef core::pose::metrics::StructureDependentCalculator parent;

	/// @brief ctor for positions, dist_cutoff will be initialized using default value from option system
	NeighborhoodByDistanceCalculator(std::set< core::Size > const & central_residues);

	/// @brief ctor for positions, with custom dist_cutoff supplied by user
	NeighborhoodByDistanceCalculator(std::set< core::Size > const & central_residues, core::Real dist_cutoff);

	/// @brief copy ctor
	NeighborhoodByDistanceCalculator( NeighborhoodByDistanceCalculator const & calculator );

	virtual core::pose::metrics::PoseMetricCalculatorOP clone() const;

	//accessors for constant/input data
	/// @brief return central residues set
	std::set< core::Size > const & central_residues() const { return central_residues_; }

	/// @brief return distance cutoff
	core::Real dist_cutoff() const { return dist_cutoff_; }

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & pose );

private:

	/// @brief whose neighbors are we finding?
	std::set< core::Size > central_residues_;
	/// @brief stores the input - how far away is a neighbor?
	core::Real const dist_cutoff_;
	/// @brief the number of neighbors, INCLUSIVE of central residues
	core::Size num_neighbors_;
	/// @brief the number of neighbors for each of the central residues
	std::map<core::Size, core::Size> num_neighbors_map_;
	/// @brief the set of neighbors, INCLUSIVE of central_residues
	std::set< core::Size > neighbors_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NeighborhoodByDistanceCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_NeighborhoodByDistanceCalculator )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_toolbox_PoseMetricCalculators_NeighborhoodByDistanceCalculator_HH

