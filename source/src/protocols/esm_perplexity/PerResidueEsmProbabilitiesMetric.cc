// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/esm_perplexity/PerResidueEsmProbabilitiesMetric.cc
/// @brief A class for predicting amino acid probabilities P(Sequence|Sequence) using the ESM language model.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

// Unit headers
#include <protocols/esm_perplexity/PerResidueEsmProbabilitiesMetric.hh>
#include <protocols/esm_perplexity/PerResidueEsmProbabilitiesMetricCreator.hh>
#include <protocols/esm_perplexity/EsmPerplexityTensorflowProtocol.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>


// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.esm_perplexity.PerResidueEsmProbabilitiesMetric" );


namespace protocols {
namespace esm_perplexity {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueEsmProbabilitiesMetric::PerResidueEsmProbabilitiesMetric():
	core::simple_metrics::PerResidueProbabilitiesMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueEsmProbabilitiesMetric::~PerResidueEsmProbabilitiesMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueEsmProbabilitiesMetric::PerResidueEsmProbabilitiesMetric( PerResidueEsmProbabilitiesMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueEsmProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< PerResidueEsmProbabilitiesMetric >( *this );
}

std::string
PerResidueEsmProbabilitiesMetric::name() const {
	return name_static();
}

std::string
PerResidueEsmProbabilitiesMetric::name_static() {
	return "PerResidueEsmProbabilitiesMetric";

}
std::string
PerResidueEsmProbabilitiesMetric::metric() const {
	return name_static();
}

void
PerResidueEsmProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	model_ = tag->getOption< std::string >( "model" );

	write_pssm_ = tag->getOption< std::string >( "write_pssm", "" );

	multirun_ = tag->getOption< bool >( "multirun", true );
#ifdef USE_TENSORFLOW
	// checking that we have fitting floating point values
	runtime_assert_string_msg( sizeof( float ) == 4 && CHAR_BIT == 8, "Model depends on 32-bit floating point values, if you see this message your operating system doesn't support this" );
	tensorflow_protocol_ = utility::pointer::make_shared< EsmPerplexityTensorflowProtocol > ( model_ );
#else
	utility_exit_with_message( "The PerResidueEsmProbabilitiesMetric needs to be compiled with either extras=tensorflow or extras=tensorflow_gpu! ");
#endif // USE_TENSORFLOW
	set_residue_selector(
		core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" )
	);
	if ( tag->hasOption("attention_mask_selection") ) {
		set_residue_selector2( core::select::residue_selector::parse_residue_selector( tag, datamap, "attention_mask_selection"));
	}
}

void
PerResidueEsmProbabilitiesMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "model", xs_string, "Which ESM model to use for the prediction" );

	attlist + XMLSchemaAttribute::attribute_w_default( "write_pssm", xs_string, "Output filename for the psi-blast like position-specific-scoring-matrix to be used with the FavorSequenceProfile Mover", "" );

	attlist + XMLSchemaAttribute::attribute_w_default( "multirun", xsct_rosetta_bool, "Whether to run multirun the network (one inference pass for all selected residues", "true");

	attributes_for_parse_residue_selector( attlist, "residue_selector", "A residue selector specifying which residue or residues to predict on." );
	attributes_for_parse_residue_selector( attlist, "attention_mask_selection", "A residue selector specifying which residues to mask.");


	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for estimating the probability of an amino acid at a given position, as predicted by the ESM language model.", attlist);
}

std::map< core::Size, std::map< core::chemical::AA, core::Real >>
PerResidueEsmProbabilitiesMetric::calculate( core::pose::Pose const &pose ) const {

	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );
	if ( residue_selector_ != nullptr ) {
		selection = residue_selector_->apply(pose);
	}
	core::select::residue_selector::ResidueSubset attention_mask_selection( pose.total_residue(), false );
	if ( selector_two_ != nullptr ) {
		attention_mask_selection = selector_two_->apply(pose);
	}
	std::map< core::Size, std::map< core::chemical::AA, core::Real >> return_probabilities_map;
#ifdef USE_TENSORFLOW
    // predict the logits
    std::map< core::Size, utility::vector1< core::Real >> logit_map = tensorflow_protocol_->compute_logits(pose, selection, attention_mask_selection, multirun_ );
    // write a PSSM in psi-blast format, we need to fill a map with logits instead of probabilities first
    if ( !write_pssm_.empty() ) {
        // first get the sequence of all positions that were predicted on
        std::string selection_sequence = get_selection_sequence( pose.sequence(), logit_map );
        std::map< core::Size, std::map< core::chemical::AA, core::Real >> logit_aa_map;
        fill_return_map( logit_map, logit_aa_map );
        output_sequence_profile( selection_sequence, logit_aa_map, write_pssm_ );
    }
    // convert to probabilities
    std::map< core::Size, utility::vector1< core::Real >> softmax_map;
    EsmPerplexityTensorflowProtocol::softmax( logit_map, softmax_map );
    // only get the probabilities for amino acids excluding special tokens
    fill_return_map( softmax_map, return_probabilities_map );
#else
	utility_exit_with_message( "This metric needs to be compiled with either extras=tensorflow or extras=tensorflow_gpu! " );
#endif // USE_TENSORFLOW
	return return_probabilities_map;
}

/// @brief Fill the return_map with the probabilities from the softmax_map for all amino acids
/// @param[in] softmax_map The map containing all probabilities or logits returned from ESM (including special tokens)
/// @param[in] return_map The map that will be returned by the calculate function, will be filled with probabilities by this function
void
PerResidueEsmProbabilitiesMetric::fill_return_map( std::map< core::Size, utility::vector1< core::Real >> const &softmax_map, std::map< core::Size, std::map< core::chemical::AA, core::Real>> &return_map ) {
	// get the probabilities for all amino acids
	for ( const auto &pair: softmax_map ) {
		core::Size selected_position = pair.first;
		utility::vector1<core::Real> softmax_vec = pair.second;
		// we skip the first 4 values as those are special characters
		for ( core::Size prob_index = 4;
				prob_index < EsmPerplexityTensorflowProtocol::alphabet_.size(); ++prob_index ) {
			char const curr_aa = EsmPerplexityTensorflowProtocol::alphabet_[prob_index];
			core::chemical::AA const aa_enum = core::chemical::aa_from_oneletter_code(curr_aa);
			core::Real probability = softmax_vec[prob_index + 1];
			// check for NaN/inf to avoid problems later
			if ( std::isnan(probability) || std::isinf(probability) ) {
				probability = 0.0;
			}

			return_map[selected_position][aa_enum] = probability;
		}
	}
}

/// @brief Get the sequence of the selection
std::string
PerResidueEsmProbabilitiesMetric::get_selection_sequence( std::string const & pose_sequence, std::map< core::Size, utility::vector1< core::Real >> const & position_map) const {
	std::string selection_sequence;
	for ( auto const & pair : position_map ) {
		selection_sequence += pose_sequence[pair.first - 1];
	}
	return selection_sequence;
}


/// @brief Set the residue selector that we'll be using.
/// @details Passing nullptr results in no residue selector being used.
void
PerResidueEsmProbabilitiesMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = std::move(selector_in);
}
/// @brief set the optional residue selector for the attention masking
void
PerResidueEsmProbabilitiesMetric::set_residue_selector2(core::select::residue_selector::ResidueSelectorCOP selector) {
	selector_two_ = std::move(selector);
}
/// @brief Get the residue selector.
/// @details If this returns nullptr, it means that no residue selector is being used.
core::select::residue_selector::ResidueSelectorCOP
PerResidueEsmProbabilitiesMetric::residue_selector() const {
	return residue_selector_;
}

void
PerResidueEsmProbabilitiesMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueEsmProbabilitiesMetric::provide_xml_schema( xsd );
}

std::string
PerResidueEsmProbabilitiesMetricCreator::keyname() const {
	return PerResidueEsmProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueEsmProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PerResidueEsmProbabilitiesMetric >();
}

void
PerResidueEsmProbabilitiesMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PerResidueEsmProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the PerResidueEsmProbabilitiesMetric."
		)
	);
	citations.add( EsmPerplexityTensorflowProtocol::get_ESM_neural_net_citation() );
}
} //esm_perplexity
} //protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::esm_perplexity::PerResidueEsmProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric>( this ) );
	arc( CEREAL_NVP ( model_ ) );
	arc( CEREAL_NVP ( write_pssm_ ) );
	arc( CEREAL_NVP ( multirun_ ) );
	arc( CEREAL_NVP (residue_selector_) );
	arc( CEREAL_NVP (selector_two_));
	// EXEMPT tensorflow_protocol_
}

template< class Archive >
void
protocols::esm_perplexity::PerResidueEsmProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric >( this ) );
	arc( model_ );
	arc( write_pssm_ );
	arc( multirun_ );
	core::select::residue_selector::ResidueSelectorOP residue_selector_local_;
	arc( residue_selector_local_);
	residue_selector_ = residue_selector_local_;

	tensorflow_protocol_ = utility::pointer::make_shared< EsmPerplexityTensorflowProtocol > ( model_ );

	std::shared_ptr< core::select::residue_selector::ResidueSelector > local_selector;
	arc( local_selector );
	selector_two_ = local_selector;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::esm_perplexity::PerResidueEsmProbabilitiesMetric );
CEREAL_REGISTER_TYPE( protocols::esm_perplexity::PerResidueEsmProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_esm_perplexity_PerResidueEsmProbabilitiesMetric )
#endif // SERIALIZATION





