// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/MinimumInterAlphaDistanceMetric.cc
/// @brief min dist
/// @author Frank Teets (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/simple_metrics/MinimumInterAlphaDistanceMetric.hh>
#include <protocols/pose_sewing/simple_metrics/MinimumInterAlphaDistanceMetricCreator.hh>

// Core headers
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

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

#include <numeric/HomogeneousTransform.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

static basic::Tracer TR( "protocols.pose_sewing.simple_metrics.MinimumInterAlphaDistanceMetric" );


namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
MinimumInterAlphaDistanceMetric::MinimumInterAlphaDistanceMetric():
	core::simple_metrics::RealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
MinimumInterAlphaDistanceMetric::~MinimumInterAlphaDistanceMetric(){}

core::simple_metrics::SimpleMetricOP
MinimumInterAlphaDistanceMetric::clone() const {
	return utility::pointer::make_shared< MinimumInterAlphaDistanceMetric >( *this );
}

std::string
MinimumInterAlphaDistanceMetric::name() const {
	return name_static();
}

std::string
MinimumInterAlphaDistanceMetric::name_static() {
	return "MinimumInterAlphaDistanceMetric";

}
std::string
MinimumInterAlphaDistanceMetric::metric() const {

	return "SHORT_NAME_FOR_SCOREFILE_HEADER_DEFAULT";
}

void
MinimumInterAlphaDistanceMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &  datamap)
{
	SimpleMetric::parse_base_tag( tag );

	if ( tag->hasOption("sequence_gap") ) {
		sequence_gap_ = tag->getOption<core::Size>("sequence_gap");
	}

	if ( tag->hasOption( "A_selector" ) ) {
		A_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "A_selector" ), datamap );
	}

	if ( tag->hasOption( "B_selector" ) ) {
		B_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "B_selector" ), datamap );
	}

	if ( tag->hasOption( "selector" ) ) {
		A_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
		B_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}

	if ( tag->hasOption( "all_atom_selector" ) ) {
		all_atom_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "all_atom_selector" ), datamap );
	}
}

void
MinimumInterAlphaDistanceMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "sequence_gap", xsct_non_negative_integer, "minimal gap between scored residues in sequence space. Anything closer than this in primary sequence will be ignored.","6" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector",  "Residue selector for all-v-all comparison" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "A_selector",  "first residue selector for a-v-b comparison" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "B_selector",  "second residue selector for a-v-b comparison" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "all_atom_selector",  "Replaces B selector and considers all atoms, not just alpha carbons, to the alpha carbons of A" );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"min dist", attlist);
}

core::Real
MinimumInterAlphaDistanceMetric::calculate(const core::pose::Pose & pose) const {

	core::Real min_dist = 99999;
	core::Real curr_dist = 99999;

	core::select::residue_selector::ResidueSubset A_selection( pose.total_residue(), false );
	if ( A_selector_ != nullptr ) {
		A_selection = A_selector_->apply( pose );
	}
	core::select::residue_selector::ResidueSubset B_selection( pose.total_residue(), false );
	if ( B_selector_ != nullptr ) {
		B_selection = B_selector_->apply( pose );
	}
	core::select::residue_selector::ResidueSubset all_atom_selection( pose.total_residue(), false );
	if ( all_atom_selector_ != nullptr ) {
		all_atom_selection = all_atom_selector_->apply( pose );
	}
	for ( core::Size n_res = 1; n_res <= pose.size(); n_res++ ) {
		if ( A_selection[n_res] ) {
			numeric::xyzVector<core::Real> nterm_CA = pose.residue(n_res).xyz(2);
			for ( core::Size c_res = n_res+sequence_gap_; c_res <= pose.size(); c_res++ ) {
				if ( B_selection[c_res] ) {
					numeric::xyzVector<core::Real> cterm_CA = pose.residue(c_res).xyz(2);
					curr_dist = nterm_CA.distance(cterm_CA);
					if ( curr_dist < min_dist ) {
						TR << curr_dist << std::endl;
						min_dist = curr_dist;
					}
				}
				if ( all_atom_selection[c_res] ) {
					for ( core::Size curr_index = 1; curr_index <= pose.residue(c_res).natoms(); curr_index++ ) {
						numeric::xyzVector<core::Real> cterm_CA = pose.residue(c_res).xyz(curr_index);
						curr_dist = nterm_CA.distance(cterm_CA);
						if ( curr_dist < min_dist ) {
							TR << curr_dist << std::endl;
							min_dist = curr_dist;
						}
					}
				}
			}
		}
	}

	return min_dist;
}

void
MinimumInterAlphaDistanceMetric::set_selectors( core::select::residue_selector::ResidueSelectorCOP A_selector, core::select::residue_selector::ResidueSelectorCOP B_selector, core::select::residue_selector::ResidueSelectorCOP all_atom_selector){
	A_selector_ = A_selector;
	B_selector_ = B_selector;
	all_atom_selector_ = all_atom_selector;
}


/// @brief Provide authorship information for an unpublished Rosetta module.
void
MinimumInterAlphaDistanceMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"MinimumInterAlphaDistanceMetric", basic::citation_manager::CitedModuleType::Mover,
		"Frank Teets",
		"Institute for Protein Innovation",
		"frank.teets@proteininnovation.org",
		"Wrote the MinimumInterAlphaDistanceMetric."
		)
	);
}


void
MinimumInterAlphaDistanceMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	MinimumInterAlphaDistanceMetric::provide_xml_schema( xsd );
}

std::string
MinimumInterAlphaDistanceMetricCreator::keyname() const {
	return MinimumInterAlphaDistanceMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
MinimumInterAlphaDistanceMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< MinimumInterAlphaDistanceMetric >();
}

} //simple_metrics
} //pose_sewing
} //protocols




