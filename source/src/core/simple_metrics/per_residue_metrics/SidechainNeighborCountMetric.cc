// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetric.cc
/// @brief A metric for calculating each sidechains neighbors based on cones.  This metric uses the same core code as the LayerSelector.  It can be combined with the SimpleMetricSelector to select any set of residues based on burial.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetric.hh>
#include <core/simple_metrics/per_residue_metrics/SidechainNeighborCountMetricCreator.hh>
#include <core/select/util/burial_utilities.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.SidechainNeighborCountMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SidechainNeighborCountMetric::SidechainNeighborCountMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SidechainNeighborCountMetric::~SidechainNeighborCountMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SidechainNeighborCountMetric::SidechainNeighborCountMetric( SidechainNeighborCountMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
SidechainNeighborCountMetric::clone() const {
	return utility::pointer::make_shared< SidechainNeighborCountMetric >( *this );
}

std::string
SidechainNeighborCountMetric::name() const {
	return name_static();
}

std::string
SidechainNeighborCountMetric::name_static() {
	return "SidechainNeighborCountMetric";

}
std::string
SidechainNeighborCountMetric::metric() const {
	return "sc_nbr_counts";
}

/// @brief Set the exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.
/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
/// vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the distance factor.
void
SidechainNeighborCountMetric::set_angle_exponent( core::Real angle_exponent ){
	angle_exponent_ = angle_exponent;
}

/// @brief Set the shift factor for the angle term.
/// @details Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom
/// vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance factor.
void
SidechainNeighborCountMetric::set_angle_shift_factor( core::Real angle_shift_factor ){
	angle_shift_factor_ = angle_shift_factor;
}

/// @brief Set the exponent for the distance term, which affects how sharp the falloff is with distance.
/// @details The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,
/// and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.
void
SidechainNeighborCountMetric::set_dist_exponent( core::Real dist_exponent){
	dist_exponent_ = dist_exponent;
}

/// @brief Midpoint of the distance falloff sigmoid.
/// @details Defaults to 9.0.  Only used by the sidchain_neighbors code.
void
SidechainNeighborCountMetric::set_dist_midpoint( core::Real dist_midpoint ){
	dist_midpoint_ = dist_midpoint;
}

/// @brief Factor by which number of residue neighbors is divided.
/// @details Defaults to 1.0.
void
SidechainNeighborCountMetric::set_res_denominator( core::Real res_denom ){
	res_denominator_ = res_denom;
}

void
SidechainNeighborCountMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

	set_angle_exponent( tag->getOption< core::Real >("angle_exponent", angle_exponent_ ));
	set_angle_shift_factor( tag->getOption< core::Real >("angle_shift_factor", angle_shift_factor_));
	set_dist_exponent( tag->getOption< core::Real >("dist_exponent", dist_exponent_) );
	set_dist_midpoint( tag->getOption< core::Real >("dist_midpoint", dist_midpoint_) );
	set_res_denominator( tag->getOption< core::Real >("res_denominator", res_denominator_));

}

void
SidechainNeighborCountMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	std::string docs = "The exponent for the angle term, which affects how close other atoms have to be to the CA-CB line to be counted fully.\n"
		"Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom"
		" vector.  The shift factor is then added, and the resulting value is raised to the angle_exponent (default 2.0) and multiplied by the"
		" distance factor.";

	attlist + XMLSchemaAttribute::attribute_w_default("angle_exponent", xsct_real, docs, "2.0");

	docs = "The shift factor for the angle term.\n"
		"Residues in the cone are counted and the count is multiplied by the cosine of the angle between the CA-CB vector and the CA-other atom"
		" vector.  The shift factor (default 0.5) is then added, and the resulting value is raised to the angle_exponent and multiplied by the distance"
		" factor.";
	attlist + XMLSchemaAttribute::attribute_w_default("angle_shift_factor", xsct_real, docs, "0.5");

	docs = "Set the exponent for the distance term, which affects how sharp the falloff is with distance.\n"
		"The distance term is: 1/(1+exp(n*(d - m))), where d is the distance, n is the exponent set by this term,"
		" and m is the midpoint of the falloff.  The n value sets the sharpness.  Defaults to 1.0.";

	attlist + XMLSchemaAttribute::attribute_w_default("dist_exponent", xsct_real, docs, "1.0");

	docs = "Midpoint of the distance falloff sigmoid.\n"\
		" Defaults to 9.0.  Only used by the sidchain_neighbors code.";

	attlist + XMLSchemaAttribute::attribute_w_default("dist_midpoint", xsct_real, docs, "9.0");

	docs = "Factor by which number of residue neighbors is divided.\n"\
		" Defaults to 1.0.";

	attlist + XMLSchemaAttribute::attribute_w_default("res_denominator", xsct_real, docs, "1.0");



	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		" A metric for calculating each sidechains neighbors based on cones.  This metric uses the same core code as the LayerSelector.  It can be combined with the SimpleMetricSelector to select any set of residues based on burial.", attlist);
}

std::map< core::Size, core::Real >
SidechainNeighborCountMetric::calculate(const pose::Pose & pose) const {
	std::map< core::Size, core::Real > final_counts;

	utility::vector1< core::Size > counts = core::select::util::calc_sc_neighbors(pose, angle_exponent_, angle_shift_factor_, dist_exponent_, dist_midpoint_, res_denominator_);

	assert(pose.size() == counts.size());

	for ( core::Size res : core::select::get_residues_from_subset( get_selector()->apply( pose )) ) {
		final_counts[res] = counts[res];
	}
	return final_counts;
}

void
SidechainNeighborCountMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SidechainNeighborCountMetric::provide_xml_schema( xsd );
}

std::string
SidechainNeighborCountMetricCreator::keyname() const {
	return SidechainNeighborCountMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SidechainNeighborCountMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< SidechainNeighborCountMetric >();
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION

core::Real angle_exponent_ = 2.0;
core::Real angle_shift_factor_ = 0.5;
core::Real dist_exponent_ = 1.0;
core::Real dist_midpoint_ = 9.0;
core::Real res_denominator_ = 1.0;

template< class Archive >
void
core::simple_metrics::per_residue_metrics::SidechainNeighborCountMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( angle_exponent_ ) );
	arc( CEREAL_NVP( angle_shift_factor_));
	arc( CEREAL_NVP( dist_exponent_));
	arc( CEREAL_NVP( dist_midpoint_ ));
	arc( CEREAL_NVP( res_denominator_));

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::SidechainNeighborCountMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( angle_exponent_ );
	arc( angle_shift_factor_);
	arc( dist_exponent_);
	arc( dist_midpoint_);
	arc( res_denominator_);
}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::SidechainNeighborCountMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::SidechainNeighborCountMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_SidechainNeighborCountMetric )
#endif // SERIALIZATION




