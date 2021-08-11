// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/BlockwisePoseCompMotifMetric.cc
/// @brief composite metric returning motif scores that pass per block
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/simple_metrics/BlockwisePoseCompMotifMetric.hh>
#include <protocols/pose_sewing/simple_metrics/BlockwisePoseCompMotifMetricCreator.hh>
#include <protocols/pose_sewing/util.hh>

// Core headers
#include <core/simple_metrics/CompositeRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/pose/Pose.hh>

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
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <protocols/sewing/scoring/MotifScorer.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/selection.hh>
#include <protocols/pose_sewing/data_storage/DsspShiftArray.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>
#include <core/select/residue_selector/BlockSelector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.pose_sewing.simple_metrics.BlockwisePoseCompMotifMetric" );


namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
BlockwisePoseCompMotifMetric::BlockwisePoseCompMotifMetric():
	core::simple_metrics::CompositeRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BlockwisePoseCompMotifMetric::~BlockwisePoseCompMotifMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
BlockwisePoseCompMotifMetric::BlockwisePoseCompMotifMetric( BlockwisePoseCompMotifMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
BlockwisePoseCompMotifMetric::clone() const {
	return utility::pointer::make_shared< BlockwisePoseCompMotifMetric >( *this );
}

std::string
BlockwisePoseCompMotifMetric::name() const {
	return name_static();
}

std::string
BlockwisePoseCompMotifMetric::name_static() {
	return "BlockwisePoseCompMotifMetric";

}
std::string
BlockwisePoseCompMotifMetric::metric() const {
	return "BlockwisePoseCompMotifMetric";

}

utility::vector1< std::string >
BlockwisePoseCompMotifMetric::get_metric_names() const {
	utility::vector1< std::string > dummy_out;
	return dummy_out;
}

void
BlockwisePoseCompMotifMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}




	normalize_by_residues_ = tag->getOption<bool>("normalize_by_residues", normalize_by_residues_);
	drop_best_ = tag->getOption<bool>("drop_best", drop_best_);

	if ( tag->hasOption( "distance_mode" ) ) {
		distance_mode_ = tag->getOption<bool>("distance_mode", distance_mode_);
	}


	if ( (! selector_ ) ) {
		utility_exit_with_message("BlockwisePoseCompat: Please pass a selector!");
	}

}

void
BlockwisePoseCompMotifMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "normalize_by_residues", xs_boolean, "normalize by the number of residues rather than the residue count","false" )
		+ XMLSchemaAttribute::attribute_w_default( "drop_best", xs_boolean, "drop_single_best_interaction","false" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_mode", xs_boolean, "use residue distances instead of motifs","false" );


	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector",  "Selector to calculate Blocks (required) " );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"Calculate motif scores between blocks or elements. ", attlist);
}

std::map< std::string, core::Real >
BlockwisePoseCompMotifMetric::calculate(const core::pose::Pose & pose) const {

	std::set<core::select::residue_selector::ResidueSubset> block_selections;
	if ( ! selector_ ) {
		utility_exit_with_message("BlockwiseShapeComp requires a selector for blockwise calculations!");
	}
	core::select::residue_selector::ResidueSubset selection = selector_->apply( pose );
	identify_ss_blocks(block_selections, selection );

	std::map< std::string, core::Real > out;
	TR << "found " << block_selections.size() << std::endl;
	if ( block_selections.size() < 2 ) {
		TR << "Blockwise invalid with <2 blocks!" << std::endl;
		return out;
	}
	calculate_bw_pose_compat_motifs(out, pose, block_selections, drop_best_, normalize_by_residues_,-2.0,!distance_mode_);
	return out;
}

void
BlockwisePoseCompMotifMetric::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
BlockwisePoseCompMotifMetric::set_normalize_by_residues(bool in){
	normalize_by_residues_ = in;
}

void
BlockwisePoseCompMotifMetric::set_drop_best(bool in){
	drop_best_ = in;
}

void
BlockwisePoseCompMotifMetric::set_distance_mode(bool in){
	distance_mode_ = in;
}

void
BlockwisePoseCompMotifMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	BlockwisePoseCompMotifMetric::provide_xml_schema( xsd );
}

std::string
BlockwisePoseCompMotifMetricCreator::keyname() const {
	return BlockwisePoseCompMotifMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
BlockwisePoseCompMotifMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< BlockwisePoseCompMotifMetric >( );
}

} //simple_metrics
} //pose_sewing
} //protocols





