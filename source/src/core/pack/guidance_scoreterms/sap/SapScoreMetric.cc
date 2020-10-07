// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapScoreMetric.cc
/// @brief A metric to report the sap score of a selected region.
/// @detail Developability index: a rapid in silico tool for the screening of antibody aggregation propensity
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapScoreMetric.hh>
#include <core/pack/guidance_scoreterms/sap/SapScoreMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/guidance_scoreterms/sap/util.hh>

#include <core/scoring/xml_util.hh>
#include <core/pose/Pose.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.pack.guidance_scoreterms.sap.SapScoreMetric" );


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


/////////////////////
/// Constructors  ///
/////////////////////

SapScoreMetric::SapScoreMetric(
	core::select::residue_selector::ResidueSelectorCOP score_selector /*= TrueSelector*/,
	core::select::residue_selector::ResidueSelectorCOP sap_calculate_selector /*= nullptr*/,
	core::select::residue_selector::ResidueSelectorCOP sasa_selector /*= nullptr*/
) :
	core::simple_metrics::RealMetric()
{
	set_score_selector( score_selector );
	set_sap_calculate_selector( sap_calculate_selector );
	set_sasa_selector( sasa_selector );
}


SapScoreMetric::~SapScoreMetric(){}


SapScoreMetric::SapScoreMetric( SapScoreMetric const & ot ) :
	core::simple_metrics::RealMetric( ot )
{
	*this = ot;
}


SapScoreMetric &
SapScoreMetric::operator=( SapScoreMetric const & ot ) {
	core::simple_metrics::RealMetric::operator=( ot );
	runtime_assert( ot.score_selector_ );
	set_score_selector( ot.score_selector_ );
	set_sap_calculate_selector( ot.sap_calculate_selector_ );
	set_sasa_selector( ot.sasa_selector_ );
	return *this;
}



core::simple_metrics::SimpleMetricOP
SapScoreMetric::clone() const {
	return utility::pointer::make_shared<SapScoreMetric>( *this );
}

void
SapScoreMetric::set_score_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ){
	if ( ! selector ) {
		utility_exit_with_message("SapScoreMetric: You can't pass a nullptr to set_score_selector(). It defaults to TrueSelector anyways.");
	}
	score_selector_ = selector->clone();
}

void
SapScoreMetric::set_sap_calculate_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ){
	if ( ! selector ) {
		sap_calculate_selector_ = nullptr;
		return;
	}
	sap_calculate_selector_ = selector->clone();
}

void
SapScoreMetric::set_sasa_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ){
	if ( ! selector ) {
		sasa_selector_ = nullptr;
		return;
	}
	sasa_selector_ = selector->clone();
}


std::string
SapScoreMetric::name() const {
	return name_static();
}

std::string
SapScoreMetric::name_static() {
	return "SapScoreMetric";

}
std::string
SapScoreMetric::metric() const {
	return "sap_score";
}

void
SapScoreMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption("score_selector") ) {
		set_score_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "score_selector" ) );
	}
	if ( tag->hasOption("sap_calculate_selector") ) {
		set_sap_calculate_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "sap_calculate_selector" ) );
	}
	if ( tag->hasOption("sasa_selector") ) {
		set_sasa_selector( core::select::residue_selector::parse_residue_selector( tag, datamap, "sasa_selector" ) );
	}
}

void
SapScoreMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"score_selector", xs_string,
		"Which residues should be included in the sap score? Optional, will default to full-pose.",
		"true_selector" )
		+ XMLSchemaAttribute(
		"sap_calculate_selector", xs_string,
		"Which residues should be present during the sap calculation? Only residues in the score_selector will have their values reported, but residues"
		" in this selector will still be assigned atom-saps which can affect the residues in score_selector. Optional, will default to score_selector." )
		+ XMLSchemaAttribute(
		"sasa_selector", xs_string,
		"Which residues should be present during the sasa calculation? Optional, will default to sap_calculate_selector." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric to report the SapScore of a specific region of a pose. Also see the AddSapConstraintMover. See this paper for more info on sap:"
		" Developability index: a rapid in silico tool for the screening of antibody aggregation propensity. Lauer, et. al. J Pharm Sci 2012",
		attlist);
}

core::Real
SapScoreMetric::calculate( core::pose::Pose const & pose ) const {
	runtime_assert( score_selector_ );

	core::select::residue_selector::ResidueSelectorCOP sap_calculate_sel = bool(sap_calculate_selector_) ? sap_calculate_selector_ : score_selector_;
	core::select::residue_selector::ResidueSelectorCOP sasa_sel = bool(sasa_selector_) ? sasa_selector_ : sap_calculate_sel;

	Real sap = core::pack::guidance_scoreterms::sap::calculate_sap( pose, score_selector_, sap_calculate_sel, sasa_sel );

	return sap;
}

void
SapScoreMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SapScoreMetric::provide_xml_schema( xsd );
}

std::string
SapScoreMetricCreator::keyname() const {
	return SapScoreMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
SapScoreMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared<SapScoreMetric>();

}

} //sap
} //guidance_scoreterms
} //pack
} //core



#ifdef    SERIALIZATION
template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapScoreMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );
	arc( CEREAL_NVP( score_selector_ ) );
	arc( CEREAL_NVP( sap_calculate_selector_ ) );
	arc( CEREAL_NVP( sasa_selector_ ) );
}

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapScoreMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::RealMetric >( this ) );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector );
	score_selector_ = local_selector;

	arc( local_selector );
	sap_calculate_selector_ = local_selector;

	arc( local_selector );
	sasa_selector_ = local_selector;
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::guidance_scoreterms::sap::SapScoreMetric );
CEREAL_REGISTER_TYPE( core::pack::guidance_scoreterms::sap::SapScoreMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapScoreMetric )
#endif // SERIALIZATION





