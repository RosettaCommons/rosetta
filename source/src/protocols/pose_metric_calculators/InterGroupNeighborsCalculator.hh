// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/InterGroupNeighborsCalculator.hh
/// @brief This is complicated, so pay attention.  This calculator is meant for finding interfaces between protein domains - like protein-protein interfaces but within a protein.  It's more flexible than that, though.  You define groups of residues within a protein (say, the N and C terminal domains).  You then define which pairs of groups you are interested in.  This calculator returns the union of the sets of residues at the interfaces between these domains/groups.  This calculator contains a superset of the functionality of some of the other calculators, but is less efficient in simple cases.  The pose does NOT have to have been scored.
/// @author Steven Lewis

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator_hh

//Unit headers
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/pose_metric_calculators/InterGroupNeighborsCalculator.fwd.hh>
#include <protocols/toolbox/CalcInterNeighborGroup.fwd.hh>


#include <core/pose/Pose.fwd.hh>
#include <basic/MetricValue.fwd.hh>

//Utility headers
#include <basic/options/option.hh>
#include <core/types.hh>
//#include <utility/vector1.hh>

//C++ headers
#include <set>
#include <utility> //pair

// option key includes
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {


/// @details This is complicated, so pay attention.  You define groups of residues within a protein (say, the N and C terminal domains).  You then define which pairs of groups you are interested in.  This calculator returns the union of the sets of residues at the interfaces between these domains/groups.  Functionally it is intended for "interface design" at the non-chainbreak interface between domains of multidomain proteins.  It contains a superset of the functionality of some of the other calculators (so I'll be obsoleting them, maybe?).  The pose does NOT have to have been scored.
/**
@li "groups" string returns the input groups; of type utility::vector1< std::pair< std::set< core::Size >, std::set< core::Size > > >  (not a calculated value)
@li "dist_cutoff" returns the input cutoff distance for neighbor finding (not a calculated value)
@li "neighbors" returns a std::set<core::Size> of the neighbors calculated between the group pairs.
@li "num_neighbors" returns the size of the neighbors set.
**/
class InterGroupNeighborsCalculator : public core::pose::metrics::StructureDependentCalculator {

public:
	typedef core::pose::metrics::StructureDependentCalculator parent;
	typedef std::set< core::Size > one_group;
	typedef std::pair< one_group, one_group > group_pair;
	typedef utility::vector1< group_pair > group_set;

	/// @brief
	InterGroupNeighborsCalculator(
		group_set const & groups,
		core::Real dist_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::inter_group_neighbors_cutoff]
	);

	InterGroupNeighborsCalculator( InterGroupNeighborsCalculator const & calculator );

	~InterGroupNeighborsCalculator();

	virtual core::pose::metrics::PoseMetricCalculatorOP clone() const;

	//accessors for non-recomputed input data
	/// @brief return groups
	group_set const & groups() const;

	/// @brief return distance cutoff
	core::Real dist_cutoff() const;

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & pose );

private:

	toolbox::CalcInterNeighborGroupOP calc_inter_group_;

	/// @brief stores the input - whose neighbors are we finding?
	group_set groups_;

	/// @brief stores the input - how far away is a neighbor?
	core::Real dist_cutoff_;

	/// @brief the number of neighbors in the set neighbors_
	core::Size num_neighbors_;

	/// @brief the set of neighbors to return - union of interfaces between groups
	std::set< core::Size > neighbors_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	InterGroupNeighborsCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace pose_metric_calculators
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_toolbox_PoseMetricCalculators_InterGroupNeighborsCalculator_HH
