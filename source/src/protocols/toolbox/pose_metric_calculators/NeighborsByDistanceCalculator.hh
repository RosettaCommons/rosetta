// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NeighborsByDistanceCalculator.hh
/// @brief NeighborsByDistanceCalculator can determine all the neighbors of a residue within a certain distance.  The pose does not have to have been scored (have a full Energies object).  It uses the PointGraph tools to find neighbors.  There is probably a much more sophisticated way to do this with existing Graph tools but I don't know what it is.
/// @author Steven Lewis

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborsByDistanceCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_NeighborsByDistanceCalculator_hh

//Unit headers
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.fwd.hh>


#include <core/pose/Pose.fwd.hh>
#include <basic/MetricValue.fwd.hh>

//Utility headers
#include <basic/options/option.hh>
#include <core/types.hh>

//C++ headers
#include <set>

// option key includes
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

/// @details this calculator determines the number and resids of residues within X angstroms of the given residue.  Its intended purpose is the backend for a TaskOperation that allows one to construct a PackerTask based on neighborhoods around a set of particular residues.  (It combines with a NeighborhoodByDistanceCalculator for that purpose).  It can return the identity of its seeded central residue and distance (just get functions) and calculate the neighbors and count of neighbors around that residue within that distance.  It uses the PointGraph class to do this; if you have a better/faster implementation please code it up and replace this one.  Note that returned data is INCLUSIVE of the central residue - it is part of the count and part of the std::set.
class NeighborsByDistanceCalculator : public core::pose::metrics::StructureDependentCalculator {

public:
	typedef core::pose::metrics::StructureDependentCalculator parent;

	/// @brief central_residue is the residue whose neighbors we find
	NeighborsByDistanceCalculator(
		core::Size central_residue,
		core::Real dist_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::neighbor_by_distance_cutoff]
	);

	NeighborsByDistanceCalculator( NeighborsByDistanceCalculator const & calculator );

	virtual core::pose::metrics::PoseMetricCalculatorOP clone() const;

	//accessors for non-recomputed input data
	/// @brief return central residue
	core::Size central_residue() const { return central_residue_; }

	/// @brief return distance cutoff
	core::Real dist_cutoff() const { return dist_cutoff_; }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & pose );

private:

	/// @brief stores the input - whose neighbors are we finding?
	core::Size const central_residue_;
	/// @brief stores the input - how far away is a neighbor?
	core::Real const dist_cutoff_;
	/// @brief the number of neighbors, INCLUSIVE of this residue
	core::Size num_neighbors_;
	/// @brief the set of neighbors, INCLUSIVE of this residue
	std::set< core::Size > neighbors_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NeighborsByDistanceCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_NeighborsByDistanceCalculator )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_toolbox_PoseMetricCalculators_NeighborsByDistanceCalculator_HH

