// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ScoreTermValueBasedSelector.cc
/// @brief  Selects residues based on per residue score of the given score_type.
/// Residues with scores equal to or within the lower and upper thresholds are selected.
/// @author Gerard Daniel (gerardda@uw.edu)

// Unit Headers
#include <core/select/residue_selector/ScoreTermValueBasedSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package Headers
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
//#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "core.select.residue_selector.ScoreTermValueBasedSelector" );

namespace core {
namespace select {
namespace residue_selector {

ScoreTermValueBasedSelector::ScoreTermValueBasedSelector():ResidueSelector() {}

ScoreTermValueBasedSelector::~ScoreTermValueBasedSelector() = default;

ResidueSelectorOP
ScoreTermValueBasedSelector::clone() const {
	return ResidueSelectorOP( utility::pointer::dynamic_pointer_cast<ResidueSelector>( ScoreTermValueBasedSelectorOP( new ScoreTermValueBasedSelector(*this) ) ) );
}

ResidueSubset ScoreTermValueBasedSelector::apply(const core::pose::Pose &pose) const
{
	ResidueSubset input_residues_subset( pose.total_residue(), false );
	if ( input_residues_selector_ ) {
		input_residues_subset = input_residues_selector_->apply( pose );
	} else if ( residue_nums_string_ == "ALL" ) {
		for ( core::Size i = 1; i <= input_residues_subset.size(); i++ ) {
			input_residues_subset[i] = true;
		}
	} else if ( residue_nums_string_ != "" ) {
		std::set< Size > const pose_residue_indices( get_resnum_list( residue_nums_string_, pose ) );
		for ( unsigned long pose_residue_indice : pose_residue_indices ) {
			input_residues_subset[pose_residue_indice] = true;
		}
	}

	core::pose::Pose _pose = pose;
	( *score_fxn_ )( _pose );
	ResidueSubset output_residues_subset( _pose.total_residue(), false );

	for ( core::Size res_index = 1; res_index <= _pose.total_residue(); res_index++ ) {
		if ( input_residues_subset[res_index] ) {
			core::Real weight = 0.0, score = 0.0, weighted_score = 0.0;
			// handle per residue total score here
			if ( score_type_ == core::scoring::ScoreType::end_of_score_type_enumeration ) {
				for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
					auto score_type = core::scoring::ScoreType(i);
					weight = score_fxn_->get_weight(score_type);
					if ( weight != 0 ) {
						weighted_score += weight * pose.energies().residue_total_energies(res_index)[ score_type ];
					}
				}
			} else {
				weight = score_fxn_->get_weight(score_type_);
				score = _pose.energies().residue_total_energies(res_index)[score_type_];
				weighted_score = weight * score;
			}
			if ( weighted_score >= lower_threshold_ && weighted_score <= upper_threshold_ ) {
				output_residues_subset[res_index] = true;
				TR << "Residue " << res_index << " " << pose.residue(res_index).name3() << " scores: " << weighted_score
					<< ". Selected." << std::endl;
			} else {
				TR << "Residue " << res_index << " " << pose.residue(res_index).name3() << " scores: " << weighted_score
					<< ". Not selected." << std::endl;
			}
		}
	}
	return output_residues_subset;
}

void ScoreTermValueBasedSelector::parse_my_tag(utility::tag::TagCOP tag, basic::datacache::DataMap &datamap)
{
	if ( tag->hasOption( "score_fxn" ) && datamap.has( "scorefxns", tag->getOption< std::string >( "score_fxn" ) ) ) {
		score_fxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", tag->getOption< std::string >( "score_fxn" ) );
	} else {
		score_fxn_ = core::scoring::get_score_function();
		TR << "Using default score function " << score_fxn_->get_name() << ".\n";
	}

	if ( tag->hasOption( "score_type" ) ) {
		std::string score_type = tag->getOption< std::string >( "score_type" );
		// per residue total score, special case.
		if ( score_type == "total" ) {
			score_type_ = core::scoring::ScoreType::end_of_score_type_enumeration;
		} else {
			score_type_ = core::scoring::score_type_from_name( tag->getOption< std::string >( "score_type" ) );
		}
		TR << "Using '" << score_type << "' for defining selection.\n";
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Missing score_type name! A predeclared filter should be specified as an option." );
	}

	if ( tag->hasOption( "selector" ) ) {
		input_residues_selector_ = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", tag->getOption< std::string >( "selector" ) );
		TR << "Will operate on the subset of residues resulting from selector '"
			<< input_residues_selector_->get_name() << "'." << std::endl;
	} else if ( tag->size() > 1 ) {
		//Check for embedded ResidueSelector as focus.
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "ScoreTermValueBased Residue Selector can take only a single sub residue selector tag!\n" );
		}
		input_residues_selector_ = core::select::residue_selector::get_embedded_residue_selector( tag, datamap );
		TR << "Will operate on the subset of residues resulting from "
			<< input_residues_selector_->get_name() << "." << std::endl;
	} else {
		residue_nums_string_ = tag->getOption<std::string>( "resnums", "ALL" );
		if ( residue_nums_string_ == "ALL" ) {
			TR << "Will operate on all residues." << std::endl;
		} else {
			TR << "Will operate on residues:" << residue_nums_string_ << std::endl;
		}
	}
	if ( tag->hasOption( "lower_threshold" ) ) {
		lower_threshold_ = tag->getOption< core::Real >( "lower_threshold" );
		TR << "Setting lower threshold to " << lower_threshold_
			<< ". Residues scoring equal to or above this threshold will be selected." << std::endl;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Missing lower threshold value! ScoreTermValueBasedSelector requires a lower threshold to be specified." );
	}
	if ( tag->hasOption( "upper_threshold" ) ) {
		upper_threshold_ = tag->getOption< core::Real >( "upper_threshold" );
		TR << "Setting upper threshold to " << upper_threshold_
			<< ". Residues scoring equal to or below this threshold will be selected." << std::endl;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Missing upper thershold value! ScoreTermValueBasedSelector requires a upper threshold to be specified." );
	}
}

std::string ScoreTermValueBasedSelector::get_name() const {
	return ScoreTermValueBasedSelector::class_name();
}

std::string ScoreTermValueBasedSelector::class_name() {
	return "ScoreTermValueBased";
}

void ScoreTermValueBasedSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "selector", xs_string, "residue selector to be used, in case only a subset of residues need to be considered")
		+ XMLSchemaAttribute::attribute_w_default( "resnums", xs_string, "pose number of the subset of residues to be considered. The default is 'ALL'.", "ALL" )
		+ XMLSchemaAttribute::required_attribute( "lower_threshold", xsct_real, "residues scoring equal to or above this threshold will be selected" )
		+ XMLSchemaAttribute::required_attribute( "upper_threshold", xsct_real, "residues scoring equal to or below this threshold will be selected" )
		+ XMLSchemaAttribute::required_attribute("score_type", xs_string, "the score type to be used for selection" )
		+ XMLSchemaAttribute::required_attribute("score_fxn", xs_string, "the score function to be used for evaluation" );
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), "Selects residues based on per residue score of the given score_type.", attributes );
}

ResidueSelectorOP
ScoreTermValueBasedSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new ScoreTermValueBasedSelector );
}

std::string
ScoreTermValueBasedSelectorCreator::keyname() const {
	return ScoreTermValueBasedSelector::class_name();
}

void
ScoreTermValueBasedSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ScoreTermValueBasedSelector::provide_xml_schema( xsd );
}

}
}
}
