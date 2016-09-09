// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.cc
/// @brief
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

// Utility headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/options/option.hh>

#include <utility/assert.hh>

// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

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

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {


void InterfaceNeighborDefinitionCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "list_interface" ) {
		basic::check_cast( valptr, &list_interface_, "list_interface expects to return a std::pair< utility::vector1<Size>, utility::vector1<Size> >" );
		(static_cast<basic::MetricValue<std::pair< utility::vector1<Size>, utility::vector1<Size> > > *>(valptr))->set( list_interface_ );

	} else if ( key == "interface_residues" ) {
		basic::check_cast( valptr, &interface_residues_, "interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( interface_residues_ );

	} else if ( key == "first_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain1_interface_residues_, "first_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain1_interface_residues_ );

	} else if ( key == "second_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain2_interface_residues_, "second_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain2_interface_residues_ );

	} else if ( key == "num_interface_residues" ) {
		basic::check_cast( valptr, &num_interface_residues_, "num_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_interface_residues_ );

	} else if ( key == "num_first_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain1_interface_residues_, "num_first_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain1_interface_residues_ );

	} else if ( key == "num_second_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain2_interface_residues_, "num_second_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain2_interface_residues_ );

	} else if ( key == "first_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch1_begin_num_, "first_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_begin_num_ );

	} else if ( key == "first_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch1_end_num_, "first_chain_last_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_end_num_ );

	} else if ( key == "second_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch2_begin_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_begin_num_ );

	} else if ( key == "second_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch2_end_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_end_num_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string InterfaceNeighborDefinitionCalculator::print( std::string const & key ) const {

	if ( key == "list_interface" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "first_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "second_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "num_interface_residues" ) {
		return utility::to_string( num_interface_residues_ );
	} else if ( key == "num_first_chain_interface_residues" ) {
		return utility::to_string( num_chain1_interface_residues_ );
	} else if ( key == "num_second_chain_interface_residues" ) {
		return utility::to_string( num_chain2_interface_residues_ );
	} else if ( key == "first_chain_first_resnum" ) {
		return utility::to_string( ch1_begin_num_ );
	} else if ( key == "first_chain_last_resnum" ) {
		return utility::to_string( ch1_end_num_ );
	} else if ( key == "second_chain_first_resnum" ) {
		return utility::to_string( ch2_begin_num_ );
	} else if ( key == "second_chain_last_resnum" ) {
		return utility::to_string( ch2_end_num_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


void InterfaceNeighborDefinitionCalculator::recompute( Pose const & this_pose ) {

	verify_chain_setup( this_pose );

	Length const distcut(basic::options::option[basic::options::OptionKeys::pose_metrics::interface_cutoff].value());

	interface_residues_.clear();
	chain1_interface_residues_.clear();
	chain2_interface_residues_.clear();

	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( this_pose.conformation(), *pg);
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, distcut );

	// for all nodes in chain1 == for all residues in chain 1
	utility::vector1< Size > chain1_interface, chain2_interface;
	for ( Size partner1_res = ch1_begin_num_; partner1_res <= ch1_end_num_; ++partner1_res ) {
		for ( conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( partner1_res ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex( partner1_res ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {

			// get node on other edge of that node == 2nd residue index
			Size const partner2_res = edge_iter->upper_vertex();

			// if that node(residue) is in chain 2
			if ( ( partner2_res >= ch2_begin_num_ ) && (partner2_res <= ch2_end_num_ ) ) {
				interface_residues_.insert( partner1_res );
				chain1_interface_residues_.insert( partner1_res );
				interface_residues_.insert( partner2_res );
				chain2_interface_residues_.insert( partner2_res );
				chain1_interface.push_back( partner1_res ); // add partner1 residue
				chain2_interface.push_back( partner2_res ); // add partner2 residue
			} else continue;

		} // end - for all edges of node
	} // end - for all nodes in chain1
	list_interface_ = std::make_pair( chain1_interface, chain2_interface ); // populate list_interface

	num_interface_residues_ = interface_residues_.size();
	num_chain1_interface_residues_ = chain1_interface_residues_.size();
	num_chain2_interface_residues_ = chain2_interface_residues_.size();

}

} // simple_calculators
} // metrics
} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator::InterfaceNeighborDefinitionCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator >( this ) );
	arc( CEREAL_NVP( list_interface_ ) ); // std::pair<utility::vector1<core::Size>, utility::vector1<core::Size> >
	arc( CEREAL_NVP( interface_residues_ ) ); // std::set<core::Size>
	arc( CEREAL_NVP( chain1_interface_residues_ ) ); // std::set<core::Size>
	arc( CEREAL_NVP( chain2_interface_residues_ ) ); // std::set<core::Size>
	arc( CEREAL_NVP( num_interface_residues_ ) ); // core::Size
	arc( CEREAL_NVP( num_chain1_interface_residues_ ) ); // core::Size
	arc( CEREAL_NVP( num_chain2_interface_residues_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::simple_calculators::InterfaceDefinitionCalculator >( this ) );
	arc( list_interface_ ); // std::pair<utility::vector1<core::Size>, utility::vector1<core::Size> >
	arc( interface_residues_ ); // std::set<core::Size>
	arc( chain1_interface_residues_ ); // std::set<core::Size>
	arc( chain2_interface_residues_ ); // std::set<core::Size>
	arc( num_interface_residues_ ); // core::Size
	arc( num_chain1_interface_residues_ ); // core::Size
	arc( num_chain2_interface_residues_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator );
CEREAL_REGISTER_TYPE( core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_simple_calculators_InterfaceNeighborDefinitionCalculator )
#endif // SERIALIZATION
