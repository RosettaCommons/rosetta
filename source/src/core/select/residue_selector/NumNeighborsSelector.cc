// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/NumNeighborsSelector.cc
/// @brief  The NumNeighborsSelector identifies all residues that have at least X neighbors within a distance D.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/residue_selector/NumNeighborsSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/pose/Pose.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {


NumNeighborsSelector::NumNeighborsSelector() :
	count_water_( false ),
	threshold_( 17 ),
	distance_cutoff_( 10.0 )
{}
NumNeighborsSelector::~NumNeighborsSelector() {}

/// @brief Copy constructor
///
NumNeighborsSelector::NumNeighborsSelector( NumNeighborsSelector const &src) :
	count_water_( src.count_water_ ),
	threshold_( src.threshold_ ),
	distance_cutoff_( src.distance_cutoff_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP NumNeighborsSelector::clone() const { return ResidueSelectorOP( new NumNeighborsSelector(*this) ); }

NumNeighborsSelector::NumNeighborsSelector( Size threshold, Real distance_cutoff ) :
	count_water_( false ),
	threshold_( threshold ),
	distance_cutoff_( distance_cutoff )
{}

ResidueSubset
NumNeighborsSelector::apply( core::pose::Pose const & pose ) const
{
	ResidueSubset subset( pose.size(), false );

	conformation::PointGraphOP pg( new conformation::PointGraph );
	conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );
	conformation::find_neighbors( pg, distance_cutoff_ );

	if ( count_water_ ) {
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( pg->get_vertex( ii ).num_neighbors() >= threshold_ ) {
				subset[ ii ] = true;
			}
		}
	} else {
		utility::vector1< Size > non_water_neighbor_count( pose.size(), 0 );
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( pose.residue_type(ii).aa() == chemical::aa_h2o ) continue;
			for ( conformation::PointGraph::VertexClass::UpperEdgeListIter
					iter = pg->get_vertex(ii).upper_edge_list_begin(),
					iter_end = pg->get_vertex(ii).upper_edge_list_end();
					iter != iter_end; ++iter ) {
				if ( pose.residue_type( iter->upper_vertex() ).aa() != chemical::aa_h2o ) {
					++non_water_neighbor_count[ ii ];
					++non_water_neighbor_count[ iter->upper_vertex() ];
				}
			}
			if ( non_water_neighbor_count[ ii ] >= threshold_ ) {
				subset[ ii ] = true;
			}
		}
	}
	return subset;
}

void NumNeighborsSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	count_water( tag->getOption< bool >( "count_water", false ));
	threshold( tag->getOption< Size >( "threshold", 17 ) );
	distance_cutoff_ = tag->getOption< Real >( "distance_cutoff", 10.0 );
}

std::string NumNeighborsSelector::get_name() const {
	return NumNeighborsSelector::class_name();
}

std::string NumNeighborsSelector::class_name() {
	return "NumNeighbors";
}

void
NumNeighborsSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "count_water", xsct_rosetta_bool, "XRW TO DO",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "threshold", xsct_non_negative_integer, "XRW TO DO",  "17"     )
		+ XMLSchemaAttribute::attribute_w_default(  "distance_cutoff", xsct_real, "XRW TO DO",  "10.0"   );
	xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attributes );
}


bool NumNeighborsSelector::count_water() const { return count_water_; }
Size NumNeighborsSelector::threshold() const { return threshold_; }
Real NumNeighborsSelector::distance_cutoff() const { return distance_cutoff_; }
void NumNeighborsSelector::count_water( bool setting ) { count_water_ = setting; }
void NumNeighborsSelector::threshold( Size setting ) { threshold_ = setting; }
void NumNeighborsSelector::distance_cutoff( Size setting ) { distance_cutoff_ = setting; }

ResidueSelectorOP
NumNeighborsSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new NumNeighborsSelector );
}

std::string
NumNeighborsSelectorCreator::keyname() const {
	return NumNeighborsSelector::class_name();
}

void
NumNeighborsSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NumNeighborsSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::NumNeighborsSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( count_water_ ) ); // _Bool
	arc( CEREAL_NVP( threshold_ ) ); // Size
	arc( CEREAL_NVP( distance_cutoff_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::NumNeighborsSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( count_water_ ); // _Bool
	arc( threshold_ ); // Size
	arc( distance_cutoff_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::NumNeighborsSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::NumNeighborsSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_NumNeighborsSelector )
#endif // SERIALIZATION
