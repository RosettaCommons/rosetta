// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/simple_metrics/WindowPoseCompMotifMetric.cc
/// @brief composite metric returning motif scores per window
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <protocols/pose_sewing/simple_metrics/WindowPoseCompMotifMetric.hh>
#include <protocols/pose_sewing/simple_metrics/WindowPoseCompMotifMetricCreator.hh>
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

#include <protocols/sewing/scoring/MotifScorer.hh>
#include <core/conformation/Atom.hh>
#include <core/pose/selection.hh>
#include <protocols/pose_sewing/data_storage/DsspShiftArray.hh>

#include <core/chemical/ChemicalManager.hh>

#include <core/pose/motif/reference_frames.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/motif/util.hh>
#include <core/select/residue_selector/BlockSelector.hh>

#include <algorithm>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.pose_sewing.simple_metrics.WindowPoseCompMotifMetric" );


namespace protocols {
namespace pose_sewing {
namespace simple_metrics {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
WindowPoseCompMotifMetric::WindowPoseCompMotifMetric():
	core::simple_metrics::CompositeRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
WindowPoseCompMotifMetric::~WindowPoseCompMotifMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
WindowPoseCompMotifMetric::WindowPoseCompMotifMetric( WindowPoseCompMotifMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
WindowPoseCompMotifMetric::clone() const {
	return utility::pointer::make_shared< WindowPoseCompMotifMetric >( *this );
}

std::string
WindowPoseCompMotifMetric::name() const {
	return name_static();
}

std::string
WindowPoseCompMotifMetric::name_static() {
	return "WindowPoseCompMotifMetric";

}
std::string
WindowPoseCompMotifMetric::metric() const {
	return "WindowPoseCompMotifMetric";

}

utility::vector1< std::string >
WindowPoseCompMotifMetric::get_metric_names() const {
	utility::vector1< std::string > out;
	return out;
}

void
WindowPoseCompMotifMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	window_width_ = tag->getOption<core::Size>("window_width", window_width_);
	normalize_by_residues_ = tag->getOption<bool>("normalize_by_residues", normalize_by_residues_);

	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	}
	if ( (! selector_ ) ) {
		utility_exit_with_message("WindowPoseCompat: Selector required to calculate blocks!");
	}
}

void
WindowPoseCompMotifMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "window_width", xsct_non_negative_integer, "score per residue cutoff","3" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_mode", xs_boolean, "use residue distances instead of motifs","false" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector",  "Selector to create blocks. Typically E and H." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring motif scores between blocks.  Uses a contiguous sliding window to calculate.  Returns the worst score for each window.  Used as a filter.  A score cutoff of -2 is strong, -1.5 is mid, while -.5 is typically weak.", attlist);
}

std::map< std::string, core::Real >
WindowPoseCompMotifMetric::calculate(const core::pose::Pose & pose) const {
	std::map< std::string, core::Real > out;
	std::set<core::select::residue_selector::ResidueSubset> block_selections;
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );

	if ( selector_ != nullptr ) {
		selection = selector_->apply( pose );
		identify_ss_blocks( block_selections, selection );
	} else {
		utility_exit_with_message("WindowPoseCompat requires a selector for blockwise calculations!");
	}

	calculate_bw_window_motifs(out,pose, block_selections, window_width_,!distance_mode_);
	return out;

}
void
WindowPoseCompMotifMetric::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
WindowPoseCompMotifMetric::set_window_width(core::Size in){
	window_width_ = in;
}

void
WindowPoseCompMotifMetric::set_distance_mode(bool in){
	distance_mode_ = in;
}

void
WindowPoseCompMotifMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	WindowPoseCompMotifMetric::provide_xml_schema( xsd );
}

std::string
WindowPoseCompMotifMetricCreator::keyname() const {
	return WindowPoseCompMotifMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
WindowPoseCompMotifMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< WindowPoseCompMotifMetric >( );
}

} //simple_metrics
} //pose_sewing
} //protocols





