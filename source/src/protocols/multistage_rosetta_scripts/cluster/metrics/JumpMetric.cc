// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.cc
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetricCreator.hh>
#include <protocols/multistage_rosetta_scripts/cluster/util.hh>

#include <math.h>
#include <utility/pointer/memory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <numeric/HomogeneousTransform.hh>

#include <core/pose/Pose.hh>
#include <core/kinematics/Jump.hh>

#include <utility/tag/Tag.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

JumpMetric::JumpMetric() :
	dofs_( 6, 0 )
{}

JumpMetric::JumpMetric( utility::vector1< core::Real > const & dofs ) :
	dofs_( dofs )
{}

JumpMetric::~JumpMetric(){}

platform::Real
JumpMetric::distance( JumpMetric const & other ) const {
	core::Real const angle_coeff = 1.0; //TODO

	core::Real sum_of_squares = 0;
	sum_of_squares += pow( dofs_[ 1 /*TRANS_X*/ ] - other.dofs_[ 1 ], 2 );
	sum_of_squares += pow( dofs_[ 2 /*TRANS_Y*/ ] - other.dofs_[ 2 ], 2 );
	sum_of_squares += pow( dofs_[ 3 /*TRANS_Z*/ ] - other.dofs_[ 3 ], 2 );
	sum_of_squares += angle_coeff * pow( dofs_[ 4 /*ROT_X*/ ] - other.dofs_[ 4 ], 2 );
	sum_of_squares += angle_coeff * pow( dofs_[ 5 /*ROT_Y*/ ] - other.dofs_[ 5 ], 2 );
	sum_of_squares += angle_coeff * pow( dofs_[ 6 /*ROT_Z*/ ] - other.dofs_[ 6 ], 2 );

	return sqrt( sum_of_squares );
}

void
JumpMetric::analyze( core::pose::Pose const & pose, unsigned int jumpno ){
	core::kinematics::Jump const & jump = pose.jump( jumpno );
	numeric::HomogeneousTransform< core::Real > ht( jump.get_rotation(), jump.get_translation() );
	dofs_[ 1 ] = ht.px(); //TRANS_X
	dofs_[ 2 ] = ht.py(); //TRANS_Y
	dofs_[ 3 ] = ht.pz(); //TRANS_Z

	numeric::xyzVector< core::Real > euler = ht.euler_angles_rad();
	dofs_[ 4 ] = euler.x(); //ROT_X
	dofs_[ 5 ] = euler.y(); //ROT_Y
	dofs_[ 6 ] = euler.z(); //ROT_Z
}

void JumpMetric::parse_my_tag(
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {
	unsigned int jumpno = tag->getOption< unsigned int >( "jump_number" );
	analyze( pose, jumpno );
}


void
JumpMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"jump_number",
		xsct_non_negative_integer,
		"Defines the ID of the jump being used",
		"1");

	xsd_type_definition_w_attributes( xsd, "Jump", "Stores the 6-D jump info.", attlist );
}


ClusterMetricOP
JumpMetricCreator::create_metric() const {
	return utility::pointer::make_shared< JumpMetric >();
}

void
JumpMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	return JumpMetric::provide_xml_schema( xsd );
}


} // metrics
} // cluster
} // multistage_...
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::metrics::JumpMetric::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::multistage_rosetta_scripts::cluster::ClusterMetric >( this ) );
	arc( CEREAL_NVP( dofs_ ) ); // utility::vector1<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::metrics::JumpMetric::load( Archive & arc ) {
	arc( cereal::base_class< protocols::multistage_rosetta_scripts::cluster::ClusterMetric >( this ) );
	arc( dofs_ ); // utility::vector1<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::cluster::metrics::JumpMetric );
CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::cluster::metrics::JumpMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_metrics_JumpMetric )
#endif // SERIALIZATION
