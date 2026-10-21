// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueBfactorMetric.cc
/// @brief Making a b factor per residue simple metric
/// @author tydingcw (claiborne.w.tydings@vanderbilt.edu)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueBfactorMetric.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueBfactorMetricCreator.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueBfactorMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::select;
using namespace core::select::residue_selector;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueBfactorMetric::PerResidueBfactorMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueBfactorMetric::~PerResidueBfactorMetric(){}

core::simple_metrics::SimpleMetricOP
PerResidueBfactorMetric::clone() const {
	return utility::pointer::make_shared< PerResidueBfactorMetric >( *this );
}

std::string
PerResidueBfactorMetric::name() const {
	return name_static();
}

std::string
PerResidueBfactorMetric::name_static() {
	return "PerResidueBfactorMetric";

}
std::string
PerResidueBfactorMetric::metric() const {

	return "Bfact";
}

void
PerResidueBfactorMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );


	if ( tag->hasOption("atom_type") ) {
		atom_type_ = tag->getOption<std::string>("atom_type", "CA");
	}
}

void
PerResidueBfactorMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("atom_type", xs_string, "Atom of which to get the b factor value", "CA");

	//attributes_for_parse_residue_selector( attlist, "residue_selector",
	// "Selector specifying residues." );

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		"Making a b factor per residue simple metric", attlist);
}

std::map< core::Size, core::Real >
PerResidueBfactorMetric::calculate(const core::pose::Pose & pose) const {
	utility::vector1< core::Size > selection1 = selection_positions(get_selector()->apply(pose));

	std::map< core::Size, core::Real > b_fact_map;
	for ( core::Size resi_index : selection1 ) {
		//Make sure atom name is in residue before getting the b factor
		if ( pose.residue_type(resi_index).has( atom_type_ ) ) {
			b_fact_map[resi_index] = pose.pdb_info()->bfactor(resi_index,pose.residue_type(resi_index).atom_index(atom_type_));
		}
	}
	//pose.residue();
	return b_fact_map;
}

/// @brief This simple metric is unpublished.  It returns Clay Tydings as its author.
void
PerResidueBfactorMetric::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PerResidueBfactorMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Clay Tydings",
		"Vanderbilty University",
		"claiborne.w.tydings@vanderbilt.edu",
		"Wrote the PerResidueBfactorMetric."
		)
	);
}

void
PerResidueBfactorMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueBfactorMetric::provide_xml_schema( xsd );
}

std::string
PerResidueBfactorMetricCreator::keyname() const {
	return PerResidueBfactorMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueBfactorMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PerResidueBfactorMetric >();
}

} //per_residue_metrics
} //simple_metrics
} //core


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueBfactorMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( atom_type_ ) );
	//arc( CEREAL_NVP( output_as_pdb_nums_ ) );

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueBfactorMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( CEREAL_NVP( atom_type_ ) );
	//arc( output_as_pdb_nums_ );


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::PerResidueBfactorMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::PerResidueBfactorMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_PerResidueBfactorMetric )
#endif // SERIALIZATION




