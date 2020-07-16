// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/HbondMetric.cc
/// @brief A metric to report the total h-bonds of a residue, or from a set of residues to another set of residues.  Use the SummaryMetric to get total hbonds of a selection or between selections. See the WaterMediatedBridgedHBondMetric for water-mediated h-bonds.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/HbondMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/select/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.HbondMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;
using namespace core::scoring::hbonds;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
HbondMetric::HbondMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
HbondMetric::~HbondMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
HbondMetric::HbondMetric( HbondMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
HbondMetric::clone() const {
	return utility::pointer::make_shared< HbondMetric >( *this );
}

std::string
HbondMetric::name() const {
	return name_static();
}

std::string
HbondMetric::name_static() {
	return "HbondMetric";

}
std::string
HbondMetric::metric() const {
	return "hbonds";
}

void
HbondMetric::set_residue_selector2(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_two_ = selector;
}

void
HbondMetric::set_include_self(bool include_self){
	include_self_ = include_self;
}


void
HbondMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

	set_include_self(tag->getOption< bool >("include_self", include_self_));

	if ( tag->hasOption("residue_selector2") ) {
		set_residue_selector2( select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector2" ) );
	}
}

void
HbondMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attributes_for_parse_residue_selector( attlist, "residue_selector2", "Optional Selector to measure hbonds between residues in each selection, instead of ANY between selector1 and the pose.  If NO selector is given, will calculate hbonds to all [OTHER] residues.");

	attlist + XMLSchemaAttribute::attribute_w_default("include_self", xsct_rosetta_bool, "Set to include self-self hydrogen bonds: Ex: resJ - resJ", "false");

	std::string docs =
		" A metric to report the total h-bonds residues from a selection to all [OTHER] residues, or from a set of residues to another set of residues.  If No selection is given, will report ALL vs ALL.\n"
		"\n"
		" TIPS:\n"
		"  Use the SummaryMetric to get total hbonds of a selection or total number of residues having some number of hbonds. . See the WaterMediatedBridgedHBondMetric for water-mediated h-bonds.\n"
		"\n"
		"  It is recommended to use -beta (-beta_nov16 and -genpot) as your scorefunction for better detection of hbonds.\n"
		"\n"
		"  By default does not report self-self hbonds (but this is an option).\n"
		"\n"
		" AUTHORS:\n"
		"  Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"   Citation: De-Novo Glycan Modeling in Rosetta (drafting)\n";
	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		docs, attlist);
}

std::map< core::Size, core::Real >
HbondMetric::calculate(const pose::Pose & pose) const {
	core::pose::Pose local_pose = pose;

	std::map<core::Size, core::Real> hbonds;

	//Initialize map to all selection as 0
	utility::vector1< bool > const mask = get_selector()->apply( pose );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( mask[i] ) {
			hbonds[i] = 0;
		}
	}

	utility::vector1< core::Size > const sele1 = get_residues_from_subset(mask);
	utility::vector1< core::Size > sele2;

	if ( selector_two_ ) {
		utility::vector1< bool > const mask2 = selector_two_->apply(pose);
		sele2 = get_residues_from_subset(mask2);
	} else {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( (! sele1.contains(i)) || (sele1.size() == pose.size()) ) {
				sele2.push_back(i);
			}
		}
	}

	//Create HbondSet
	//NEED copy of pose here...
	HBondSet hb_set = HBondSet(local_pose, false);

	for ( core::Size const s1 : sele1 ) {
		utility::vector1< HBondCOP > s1_hbonds = hb_set.residue_hbonds(s1);
		for ( HBondCOP hb : s1_hbonds ) {
			core::Size s2 = next_hb_res(*hb, s1);
			if ( sele2.contains(s2) ) {
				if ( s2 == s1 && ! include_self_ ) {
					continue;
				}

				hbonds[s1]+=1;
			}
		}
	}
	return hbonds;

}

void
HbondMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	HbondMetric::provide_xml_schema( xsd );
}

std::string
HbondMetricCreator::keyname() const {
	return HbondMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
HbondMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< HbondMetric >();
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::per_residue_metrics::HbondMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( include_self_ ) );
	arc( CEREAL_NVP( selector_two_ ));

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::HbondMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( include_self_ );
	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector ); // ResidueSelectorCOP
	selector_two_ = local_selector;


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::HbondMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::HbondMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_HbondMetric )
#endif // SERIALIZATION




