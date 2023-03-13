// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/BFactorSelector.hh
/// @brief  A residue selector dependent on b-factor values
/// @author AmeyaHarmalkar (harmalkar.ameya24@gmail.com)

// Unit headers
#include <core/select/residue_selector/BFactorSelector.hh>
#include <core/select/residue_selector/BFactorSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/variant_util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <utility/assert.hh>

#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.BFactorSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
BFactorSelector::BFactorSelector():
	lower_res_( 1 ),
	upper_res_( 1 ),
	residue_selector_(),
	cross_chain_boundaries_( false )
{
}

BFactorSelector::BFactorSelector(
	core::Size const lower_res,
	core::Size const upper_res,
	core::select::residue_selector::ResidueSelectorCOP const selector,
	bool cross_chain_boundaries ):
	lower_res_( lower_res ),
	upper_res_( upper_res ),
	residue_selector_( selector ),
	cross_chain_boundaries_( cross_chain_boundaries  )
{
}

/// @brief Destructor.
///
BFactorSelector::~BFactorSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//BFactorSelector::BFactorSelector(BFactorSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
BFactorSelector::clone() const {
	return utility::pointer::make_shared< BFactorSelector >( *this );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
BFactorSelector::ResidueSubset
BFactorSelector::apply( const core::pose::Pose &pose) const {

	ResidueSubset input_residues_subset( pose.total_residue(), false );
	// If residue selector is provided, create a residue subset to select from
	if ( residue_selector_ ) {
		input_residues_subset = residue_selector_->apply( pose );
	} else {
		for ( core::Size i = 1; i <= input_residues_subset.size(); i++ ) {
			input_residues_subset[i] = true;
		}
	}

	// load the pose and define all false residue subset
	core::pose::Pose real_pose = pose;
	ResidueSubset output_residues_subset( real_pose.total_residue(), false );

	// For loop across all the residues
	// To do: add specificiation for contiguous residues i.e. minimum group of 3/4 should
	// be selected for useful KIC moves
	for ( core::Size res_index = 1; res_index <= real_pose.total_residue(); res_index++ ) {

		if ( input_residues_subset[res_index] ) {

			core::Real bfactor_score = real_pose.pdb_info()->bfactor( res_index, 2 );

			if ( bfactor_score >= lower_threshold_ && bfactor_score <= upper_threshold_ ) {
				output_residues_subset[res_index] = true;
				TR << "Residue " << res_index << " " << real_pose.residue(res_index).name3() << ": Selected." << std::endl;
			} else {
				TR << "Residue " << res_index << " " << real_pose.residue(res_index).name3() << ": Not selected." << std::endl;
			}
		}
	}

	// Get contiguous segments
	// If contiguous then we run through another for loop
	ResidueRanges const ranges( output_residues_subset );
	TR << "Intervals: [ " ;
	for ( auto const & range : ranges ) {
		TR << " " << range.start() << "->" << range.stop() << ", " ;
		core::Size start = range.start();
		core::Size end = range.stop();
		core::Size count = 0;
		if ( end-start+1 >= contiguous_res_ ) { // only add residues if it satisfies the minimum contiguousness requirement
			while ( ( count < lower_res_ ) && ( start > 1 ) && !( pose::is_lower_terminus( real_pose, start ) && ! cross_chain_boundaries_ ) && ( input_residues_subset[start] ) ) {
				++count;
				--start;
			}
			TR << count << " residues added to lower terminus." << std::endl;
			count = 0;
			while ( ( count < upper_res_ ) && ( end < real_pose.size() ) && ! ( pose::is_upper_terminus( real_pose, end ) && ! cross_chain_boundaries_ ) && ( input_residues_subset[end] ) ) {
				++count;
				++end;
			}
			TR << count << " residues added to upper terminus." << std::endl;
			for ( Size r=start; r<=end; r++ ) {
				if ( output_residues_subset[r] != true ) {
					output_residues_subset[r] = true;
				}
			}
		} else {
			TR << "Residue interval doesn't satisfy minimum contiguous res." << std::endl;
			for ( Size r=start; r<=end; r++ ) {
				output_residues_subset[r] = false;
			}
		}
	}
	TR << "]" << std::endl;
	return output_residues_subset;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
BFactorSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	if ( tag->hasOption( "selector" ) ) {
		residue_selector_ = get_residue_selector(tag->getOption< std::string >( "selector" ), datamap);
		TR << "Will operate on the subset of residues resulting from selector '"
			<< residue_selector_->get_name() << "'." << std::endl;
	}

	if ( tag->hasOption( "lower_bfactor" ) ) {
		// Replace with setter
		set_lower_bfactor_threshold( tag->getOption< core::Real >( "lower_bfactor" ) );
		TR << "Setting lower bfactor threshold to " << lower_threshold_
			<< ". Residues scoring equal to or above this threshold will be selected." << std::endl;
	}

	if ( tag->hasOption( "upper_bfactor" ) ) {
		// Replace with setter
		set_upper_bfactor_threshold( tag->getOption< core::Real >( "upper_bfactor" ) );
		TR << "Setting upper bfactor threshold to " << upper_threshold_
			<< ". Residues scoring equal to or below this threshold will be selected." << std::endl;
	}

	if ( tag->hasOption( "lower" ) ) {
		set_lower_residues( tag->getOption< Size >( "lower" ) );
	}

	if ( tag->hasOption( "upper" ) ) {
		set_upper_residues( tag->getOption< Size >( "upper" ) );
	}
	if ( tag->hasOption( "cross_chain_boundaries" ) ) {
		set_cross_chain_boundaries( tag->getOption< bool >( "cross_chain_boundaries" ) );
	}

	if ( tag->hasOption( "min_contiguous_res" ) ) {
		set_min_contiguous_res( tag->getOption< Size >( "min_contiguous_res" ) );
	}

}


void
BFactorSelector::set_lower_residues( core::Size const nres )
{
	lower_res_ = nres;
}

void
BFactorSelector::set_upper_residues( core::Size const nres )
{
	upper_res_ = nres;
}

void
BFactorSelector::set_min_contiguous_res( core::Size const nres )
{
	contiguous_res_ = nres;
}

void
BFactorSelector::set_lower_bfactor_threshold( core::Real const lower_thresh )
{
	lower_threshold_ = lower_thresh;
}

void
BFactorSelector::set_upper_bfactor_threshold( core::Real const upper_thresh )
{
	upper_threshold_ = upper_thresh;
}

void
BFactorSelector::set_cross_chain_boundaries( bool cross )
{
	cross_chain_boundaries_ = cross;
}


std::string BFactorSelector::get_name() const
{
	return BFactorSelector::class_name();
}

std::string BFactorSelector::class_name()
{
	return "BFactorSelector";
}

void BFactorSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	attributes_for_parse_residue_selector(attributes,  "selector");
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "lower_bfactor", xsct_real, "residues scoring equal to or above this threshold will be selected", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default(  "upper_bfactor", xsct_real, "residues scoring equal to or below this threshold will be selected", "100.0" )
		+ XMLSchemaAttribute::attribute_w_default( "lower", xsct_non_negative_integer , "lower limit for contiguous residues (N-terminus)", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "upper", xsct_non_negative_integer , "upper limit for contiguous residues (C-terminus)", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_contiguous_res", xsct_non_negative_integer ,
		"Minimum contiguous segments acceptable i.e. residues for which contiguous stretches "
		"higher than min_contiguous_res option and satisfying the threshold criteria would be "
		"selected. Say if min_contiguous_res = 3, and the BFactorSelector selects the residues "
		"[12-15,17,18-19,23-27] but with the min_contiguous_res option the final selection would be "
		"[12-15,23-27] so that only contiguous segments of residues are selected. Useful if "
		"the output of selection is passed on to Backrub or KIC movers which need pivot and "
		"cutpoints.",
		"3" )
		+ XMLSchemaAttribute( "cross_chain_boundaries", xsct_rosetta_bool, "Allow the selector to cross chain boundaries? By default false.");
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "Select residues based on bfactor values", attributes );
}

/// @brief This residue selector is unpublished.  It returns AmeyaHarmalkar as its author.
void
BFactorSelector::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BFactorSelector", basic::citation_manager::CitedModuleType::ResidueSelector,
		"AmeyaHarmalkar",
		"Johns Hopkins University",
		"aharmal1@jhu.edu",
		"Wrote the BFactorSelector."
		)
	);
}

core::select::residue_selector::ResidueSelectorOP
BFactorSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< BFactorSelector >();
}

std::string
BFactorSelectorCreator::keyname() const {
	return BFactorSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
BFactorSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	BFactorSelector::provide_xml_schema( xsd );
}

} //residue_selector
} //select
} //core

#ifdef    SERIALIZATION

// See the serialization documentation here.  There is a script you can run.
//  https://wiki.rosettacommons.org/index.php/SerializationFAQ

template< class Archive >
void
core::select::residue_selector::BFactorSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( lower_threshold_ ) ); // Real
	arc( CEREAL_NVP( upper_threshold_ ) ); // Real
	arc( CEREAL_NVP( lower_res_ ) ); // Size
	arc( CEREAL_NVP( upper_res_ ) ); // Size
	arc( CEREAL_NVP( residue_selector_ ) ); // ResidueSelectorOP
	arc( CEREAL_NVP( cross_chain_boundaries_ ) ); // Bool
	arc( CEREAL_NVP( contiguous_res_ ) ); // Size
}

template< class Archive >
void
core::select::residue_selector::BFactorSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( lower_threshold_ ); // Real
	arc( upper_threshold_ ); // Real
	arc( lower_res_ ); // Size
	arc( upper_res_ ); // Size
	arc( residue_selector_ ); // ResidueSelectorOP
	arc( cross_chain_boundaries_ ); // Bool
	arc( contiguous_res_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::BFactorSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::BFactorSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_BFactorSelector )
#endif // SERIALIZATION
