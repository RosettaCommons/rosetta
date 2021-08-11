// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/BlockwiseShapeCompMetric.cc
/// @brief composite metric returning shape complementarities that pass per block
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/pose_sewing/simple_metrics/BlockwiseShapeCompMetric.hh>
#include <protocols/pose_sewing/simple_metrics/BlockwiseShapeCompMetricCreator.hh>
#include <protocols/pose_sewing/util.hh>

// Core headers
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/simple_filters/ShapeComplementarityFilter.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <core/pose/Pose.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>

static basic::Tracer TR( "protocols.pose_sewing.simple_metrics.BlockwiseShapeCompMetric" );


namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
BlockwiseShapeCompMetric::BlockwiseShapeCompMetric():
	core::simple_metrics::CompositeRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BlockwiseShapeCompMetric::~BlockwiseShapeCompMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
BlockwiseShapeCompMetric::BlockwiseShapeCompMetric( BlockwiseShapeCompMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
BlockwiseShapeCompMetric::clone() const {
	return utility::pointer::make_shared< BlockwiseShapeCompMetric >( *this );
}

std::string
BlockwiseShapeCompMetric::name() const {
	return name_static();
}

std::string
BlockwiseShapeCompMetric::name_static() {
	return "BlockwiseShapeCompMetric";

}
std::string
BlockwiseShapeCompMetric::metric() const {
	return "bw_sc_links";
}

utility::vector1< std::string >
BlockwiseShapeCompMetric::get_metric_names() const {
	//Blocks are not known until apply time.
	utility::vector1<std::string> out;
	out.push_back("1");
	out.push_back("2");
	return out;
}

void
BlockwiseShapeCompMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );
	filtered_sc_ = tag->getOption<core::Real>( "min_sc", filtered_sc_ );
	filtered_area_ = tag->getOption<core::Real>( "min_interface", filtered_area_ );
	filtered_d_median_ = tag->getOption<core::Real>( "max_median_dist", filtered_d_median_ );
	scale_by_length_ = tag->getOption<bool>( "scale_by_length", scale_by_length_);

	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}

	if ( ! selector_ ) {
		utility_exit_with_message("BlockwiseShapeComp: Please pass a selector");
	}
}

void
BlockwiseShapeCompMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "min_sc" , xsct_real , "Bad link if the calculated sc is less than the given value." , "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "min_interface" , xsct_real , "Bad link if the calculated interface area is less than the given value. Required." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_median_dist" , xsct_real , "Bad link if the calculated median distance between the molecular surfaces is greater than the given value." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "scale_by_length" , xs_boolean , "Scale the filtered area by block length." , "false" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector",  "Selector to calculate Blocks (required)" );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring the Shape Complementarity between blocks. Reports the number of good links for each block. A link is defined as acceptable Shape Complementarity between each block and all other blocks.", attlist);
}

std::map< std::string, core::Real >
BlockwiseShapeCompMetric::calculate(const core::pose::Pose & pose) const {
	std::map< std::string, core::Real > out;

	std::set<core::select::residue_selector::ResidueSubset> block_selections;
	core::select::residue_selector::ResidueSubset selection;

	if ( ! selector_ ) {
		utility_exit_with_message("BlockwiseShapeComp requires a selector with chosen SS for blockwise calculations!");
	}
	selection = selector_->apply( pose );
	protocols::pose_sewing::identify_ss_blocks(block_selections, selection );

	TR << "found " << block_selections.size() << std::endl;

	std::set<core::select::residue_selector::ReturnResidueSubsetSelectorCOP> block_selectors;
	for ( auto selection : block_selections ) {
		block_selectors.insert(core::select::residue_selector::ReturnResidueSubsetSelectorCOP(new core::select::residue_selector::ReturnResidueSubsetSelector(selection)));
	}

	protocols::simple_filters::ShapeComplementarityFilterOP sc_filt =  protocols::simple_filters::ShapeComplementarityFilterOP( new protocols::simple_filters::ShapeComplementarityFilter);
	sc_filt->filtered_sc(filtered_sc_);
	sc_filt->filtered_area(filtered_area_);
	sc_filt->filtered_median_distance(filtered_d_median_);
	core::Size current_block = 0;
	for ( auto it_selector_1 = block_selectors.begin(); it_selector_1 != block_selectors.end(); ++it_selector_1 ) {
		current_block++;
		core::Size good_links = 0;
		for ( auto it_selector_2 = block_selectors.begin(); it_selector_2 != block_selectors.end(); ++it_selector_2 ) {
			if ( *it_selector_1 == *it_selector_2 ) {
				continue;
			}
			sc_filt->selector1(*it_selector_1);
			sc_filt->selector2(*it_selector_2);
			if ( scale_by_length_ ) {
				core::select::residue_selector::ReturnResidueSubsetSelectorCOP current_selector;
				core::Size sc_1_length = 0;
				core::Size sc_2_length = 0;
				for ( core::Size count_resnum = 1; count_resnum <= pose.size(); ++count_resnum ) {
					current_selector = *it_selector_1;
					selection=current_selector->apply(pose);
					if ( selection[count_resnum] ) {
						++sc_1_length;
					}
					current_selector = *it_selector_2;
					selection=current_selector->apply(pose);
					if ( selection[count_resnum] ) {
						++sc_2_length;
					}
				}
				if ( sc_1_length < sc_2_length ) {
					sc_filt->filtered_area(filtered_area_ * sc_1_length);
				} else {
					sc_filt->filtered_area(filtered_area_ * sc_2_length);
				}

			}
			if ( sc_filt->apply(pose) ) {
				++good_links;
			}
		}
		out[std::to_string(current_block)]= good_links;
	}
	return out;
}
void
BlockwiseShapeCompMetric::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
BlockwiseShapeCompMetric::set_filtered_sc(core::Real filtered_sc){
	filtered_sc_ = filtered_sc;
}

void
BlockwiseShapeCompMetric::set_filtered_area(core::Real filtered_area){
	filtered_area_ = filtered_area;
}

void
BlockwiseShapeCompMetric::set_filtered_d_median(core::Real filtered_d_median){
	filtered_d_median_ = filtered_d_median;
}

void
BlockwiseShapeCompMetric::set_scale_by_length(bool in){
	scale_by_length_ = in;
}

void
BlockwiseShapeCompMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	BlockwiseShapeCompMetric::provide_xml_schema( xsd );
}

std::string
BlockwiseShapeCompMetricCreator::keyname() const {
	return BlockwiseShapeCompMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
BlockwiseShapeCompMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< BlockwiseShapeCompMetric >( );
}

} //simple_metrics
} //pose_sewing
} //protocols






