// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/InterGroupNeighborsCalculator.cc
/// @brief This is complicated, so pay attention.  This calculator is meant for finding interfaces between protein domains - like protein-protein interfaces but within a protein.  It's more flexible than that, though.  You define groups of residues within a protein (say, the N and C terminal domains).  You then define which pairs of groups you are interested in.  This calculator returns the union of the sets of residues at the interfaces between these domains/groups.  This calculator contains a superset of the functionality of some of the other calculators, but is less efficient in simple cases.  The pose does NOT have to have been scored.
/// @author Steven Lewis

//Unit headers
#include <protocols/pose_metric_calculators/InterGroupNeighborsCalculator.hh>
#include <protocols/toolbox/CalcInterNeighborGroup.hh>

#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <utility/graph/Graph.hh>

//Utility headers
//#include <basic/options/option.hh>
//#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>


//C++ headers
//#include <set>
//#include <utility>

static basic::Tracer TR( "protocols.toolbox.PoseMetricCalculators.InterGroupNeighborsCalculator" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {

using namespace protocols::toolbox;

using parent = core::pose::metrics::StructureDependentCalculator;
using one_group = std::set<core::Size>;
typedef std::pair< one_group, one_group > group_pair;
using group_set = utility::vector1<group_pair>;

InterGroupNeighborsCalculator::InterGroupNeighborsCalculator( group_set const & groups, core::Real dist_cutoff )
: parent()
	//not doing anything to std::set<core::Size> - should initialize empty
{
	calc_inter_group_ = CalcInterNeighborGroupOP( new CalcInterNeighborGroup(groups, dist_cutoff) );
}

InterGroupNeighborsCalculator::InterGroupNeighborsCalculator( InterGroupNeighborsCalculator const & calculator )
: parent(), calc_inter_group_(calculator.calc_inter_group_)
{}

InterGroupNeighborsCalculator::~InterGroupNeighborsCalculator() = default;

core::pose::metrics::PoseMetricCalculatorOP InterGroupNeighborsCalculator::clone() const
{ return core::pose::metrics::PoseMetricCalculatorOP( new InterGroupNeighborsCalculator(*this) ); }

core::Real InterGroupNeighborsCalculator::dist_cutoff() const {
	return dist_cutoff_;
}

group_set const& InterGroupNeighborsCalculator::groups() const {
	return groups_;
}

void
InterGroupNeighborsCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{
	if ( key == "groups" ) {
		basic::check_cast( valptr, & groups_, "groups expects to return a utility::vector1< std::pair< std::set< core::Size >, std::set< core::Size > > >" );
		(static_cast<basic::MetricValue<group_set> *>(valptr))->set( groups_ );

	} else if ( key == "dist_cutoff" ) {
		basic::check_cast( valptr, & dist_cutoff_, "dist_cutoff expects to return a core::Real" );
		(static_cast<basic::MetricValue<core::Real> *>(valptr))->set( dist_cutoff_ );

	} else if ( key == "num_neighbors" ) {
		basic::check_cast( valptr, & num_neighbors_, "num_neighbors expects to return a core::Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( num_neighbors_ );

	} else if ( key == "neighbors" ) {
		basic::check_cast( valptr, & neighbors_, "neighbors expects to return a std::set< core::Size >" );
		(static_cast<basic::MetricValue< std::set< core::Size > > *>(valptr))->set( neighbors_ );

	} else {
		basic::Error() << "InterGroupNeighborsCalculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

} //lookup

std::string
InterGroupNeighborsCalculator::print( std::string const & key ) const
{
	if ( key == "dist_cutoff" ) {
		return utility::to_string( dist_cutoff_ );

	} else if ( key == "num_neighbors" ) {
		return utility::to_string( num_neighbors_ );

	} else if ( key == "groups" || key == "neighbors" ) {
		//set up big return string for both sets
		using namespace basic::options; //this lets you get + or (space) as spacer
		std::string const spacer( option[ OptionKeys::packing::print_pymol_selection].value() ? "+" : " ");
		std::string nbrs_string("");

		if ( key == "groups" ) {
			for ( core::Size i(1), vecsize(groups_.size()); i <= vecsize; ++i ) {
				nbrs_string += "{ (";
				for ( unsigned long it : groups_[i].first ) {
					nbrs_string += spacer + utility::to_string(it);
				}
				nbrs_string += ") ; (";
				for ( unsigned long it : groups_[i].second ) {
					nbrs_string += spacer + utility::to_string(it);
				}
				nbrs_string += ") }";
			}
			return nbrs_string;
		} else if ( key == "neighbors" ) {
			for ( unsigned long neighbor : neighbors_ ) {
				nbrs_string += utility::to_string(neighbor) + spacer;
			}
			return nbrs_string;
		}//neighbors or groups
	}//else
	basic::Error() << "InterGroupNeighborsCalculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";
} //print

void
InterGroupNeighborsCalculator::recompute( core::pose::Pose const & pose )
{
	calc_inter_group_->compute(pose);

	//JAB - Needed for check_cast
	groups_ = calc_inter_group_->groups();
	dist_cutoff_ = calc_inter_group_->dist_cutoff();
	num_neighbors_ = calc_inter_group_->num_neighbors();
	neighbors_ = calc_inter_group_->neighbors();

	return;
} //recompute

} //namespace pose_metric_calculators
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::pose_metric_calculators::InterGroupNeighborsCalculator::InterGroupNeighborsCalculator() :
	dist_cutoff_( 0 )
{}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pose_metric_calculators::InterGroupNeighborsCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( calc_inter_group_ ) );
	arc( CEREAL_NVP( groups_ ) ); // const group_set
	arc( CEREAL_NVP( dist_cutoff_ ) ); // const core::Real
	arc( CEREAL_NVP( num_neighbors_ ) ); // core::Size
	arc( CEREAL_NVP( neighbors_ ) ); // std::set<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pose_metric_calculators::InterGroupNeighborsCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( calc_inter_group_ );
	arc( const_cast< group_set & > ( groups_ ) ); // const group_set
	arc( const_cast< core::Real & > (dist_cutoff_) ); // const core::Real
	arc( num_neighbors_ ); // core::Size
	arc( neighbors_ ); // std::set<core::Size>
}
SAVE_AND_LOAD_SERIALIZABLE( protocols::pose_metric_calculators::InterGroupNeighborsCalculator );
CEREAL_REGISTER_TYPE( protocols::pose_metric_calculators::InterGroupNeighborsCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_InterGroupNeighborsCalculator )
#endif // SERIALIZATION
