// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SequenceSimilarityMetric.cc
/// @brief compares the sequences of the native structure (using native flag) to the sequence of a given pose.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <core/simple_metrics/metrics/SequenceSimilarityMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/RealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/sequence/Blosum62Map.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ref_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR( "core.simple_metrics.metrics.SequenceSimilarityMetric" );


namespace core {
namespace simple_metrics {
namespace metrics {

using namespace core::select::residue_selector;

SequenceSimilarityMetric::SequenceSimilarityMetric():
	RealMetric(),
	selector_( 0 ),
	native_pose_( 0 ),
	apply_selector_to_native_( false ),
	normalize_( true )
{}

SequenceSimilarityMetric::SequenceSimilarityMetric( select::residue_selector::ResidueSelectorCOP selector):
	RealMetric(),
	selector_( selector ),
	native_pose_( 0 ),
	apply_selector_to_native_( false ),
	normalize_( true )
{}

SequenceSimilarityMetric::~SequenceSimilarityMetric(){}

SequenceSimilarityMetric::SequenceSimilarityMetric( SequenceSimilarityMetric const & src):
	RealMetric( src )
{
	selector_ = src.selector_;
	native_pose_ = src.native_pose_;
	apply_selector_to_native_ = src.apply_selector_to_native_;
	normalize_ = src.normalize_;
}


std::string
SequenceSimilarityMetric::name() const {
	return name_static();
}

std::string
SequenceSimilarityMetric::name_static() {
	return "SequenceSimilarityMetric";

}
std::string
SequenceSimilarityMetric::metric() const {
	return "sequence_similarity";
}

core::Real
SequenceSimilarityMetric::calculate( pose::Pose const & pose ) const {

	using namespace basic::options;

	if ( native_pose_ == 0 ) {
		//std::string const filename = option[ OptionKeys::in::file::native ].value();
		//native_pose_ = import_pose::pose_from_file( filename );
		utility_exit_with_message( "Can not find native pose. If you are not using rosetta scripts, please call set_native_pose()." );
	}

	utility::vector1< bool > residues_to_count;
	if ( selector_ ) {
		if ( apply_selector_to_native_ ) {
			residues_to_count = selector_->apply( * native_pose_ );
		} else {
			residues_to_count = selector_->apply( pose );
		}
	} else {
		residues_to_count.assign( pose.size(), true );
	}

	std::string const full_pose_seq = pose.sequence();
	std::string const full_nat_seq = native_pose_->sequence();

	std::string pose_seq = "";
	std::string nat_seq  = "";

	for ( core::Size resid = 1; resid <= full_pose_seq.size() && resid <= full_nat_seq.size(); ++resid ) {
		if ( residues_to_count[ resid ] ) {
			pose_seq += full_pose_seq[ resid - 1 ];
			nat_seq += full_nat_seq[ resid - 1 ];
#ifndef NDEBUG
			//Assert that these are the same residues
			debug_assert( pose.pdb_info()->chain( resid ) == native_pose_->pdb_info()->chain( resid ) );
			debug_assert( pose.pdb_info()->number( resid ) == native_pose_->pdb_info()->number( resid ) );
#endif
		}
	}

	return score( pose_seq, nat_seq, normalize_ );
}

core::Real
SequenceSimilarityMetric::score( std::string const & seq1, std::string const & seq2, bool normalize ){
	runtime_assert( seq1.size() == seq2.size() );

	core::sequence::Blosum62Map map;

	core::Real score = 0;
	auto const num_elements = seq1.size();

	for ( Size i = 0; i < num_elements; ++i ) {
		score += map.score_for_aa_pair( seq1[ i ], seq2[ i ] );
	}

	if ( normalize ) {
		return score / num_elements;
	} else {
		return score;
	}
}

void
SequenceSimilarityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap  )
{
	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( select::residue_selector::parse_residue_selector( tag, datamap ) );
	}

	if ( datamap.has_resource("native_pose") ) {
		native_pose_ = pose::saved_native_pose( datamap );
	}

	if ( tag->hasOption( "apply_selector_to_native" ) ) {
		set_apply_selector_to_native( tag->getOption< bool >( "apply_selector_to_native" ) );
	}

	if ( tag->hasOption( "normalize" ) ) {
		set_normalize( tag->getOption< bool >( "normalize" ) );
	}

}

SimpleMetricOP
SequenceSimilarityMetric::clone() const {
	return utility::pointer::make_shared< SequenceSimilarityMetric >( *this );

}

void
SequenceSimilarityMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"apply_selector_to_native", xsct_rosetta_bool,
		"If true, apply the residue selector to the native pose instead of the given pose.",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"normalize", xsct_rosetta_bool,
		"Divide the final score by the number of positions",
		"true");

	std::string const description = "Choose which residues are to be counted." ;
	attributes_for_parse_residue_selector_default_option_name( attlist, description );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for measuring sequence similarity between the trajectry pose and the native pose. Sums the score of each selected position using BLOSUM62 (note that the score of an unmutated pose is almost certainly not 0). Must use the -native flag for this to work correctly.", attlist);
}

void
SequenceSimilarityMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SequenceSimilarityMetric::provide_xml_schema( xsd );
}

std::string
SequenceSimilarityMetricCreator::keyname() const {
	return SequenceSimilarityMetric::name_static();
}

SimpleMetricOP
SequenceSimilarityMetricCreator::create_simple_metric() const {
	return SimpleMetricOP( new SequenceSimilarityMetric );

}

} //core
} //simple_metrics
} //metrics






