// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionMetric.cc
/// @brief A class for predicting post-translational modifications using a neural nets.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adapted from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

// Unit headers
#include <protocols/ptm_prediction/PTMPredictionMetric.hh>
#include <protocols/ptm_prediction/PTMPredictionMetricCreator.hh>
#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.hh>
#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocol.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>

// Protocols headers:
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

// Citation manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.ptm_prediction.PTMPredictionMetric" );


namespace protocols {
namespace ptm_prediction {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PTMPredictionMetric::PTMPredictionMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PTMPredictionMetric::~PTMPredictionMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PTMPredictionMetric::PTMPredictionMetric( PTMPredictionMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PTMPredictionMetric::clone() const {
	return utility::pointer::make_shared< PTMPredictionMetric >( *this );
}

std::string
PTMPredictionMetric::name() const {
	return name_static();
}

std::string
PTMPredictionMetric::name_static() {
	return "PTMPredictionMetric";

}
std::string
PTMPredictionMetric::metric() const {
	return name_static();
}

void
PTMPredictionMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	modification_ = tag->getOption< std::string >( "modification" );

	ptm_mode_ = string_to_enum( modification_ );

	std::string path_to_model = "protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule/predict_PTMs/Model_" + modification_;
#ifdef USE_TENSORFLOW
// checking that we have fitting 32-bit floating point values that our tensorflow model depends on
	runtime_assert_string_msg( sizeof( float ) == 4 && CHAR_BIT == 8, "Model depends on 32-bit floating point values, if you see this message your operating system doesn't support this" );
	tensorflow_protocol_ = utility::pointer::make_shared< PTMPredictionTensorflowProtocol > ( path_to_model );
#endif // USE_TENSORFLOW
	set_residue_selector(
		core::select::residue_selector::parse_residue_selector( tag, datamap, "residue_selector" )
	);
}

void
PTMPredictionMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	std::string mod_options_string = "The post-translational modification to predict. Available options are: ";
	for ( int i = 1; i <= static_cast< int >( PTMPredictionMetricModes::N_MODES ); i++ ) {
		PTMPredictionMetricModes mod = static_cast< PTMPredictionMetricModes >( i );
		std::string mod_string = enum_to_string( mod );

		if ( i == 1 ) {
			mod_options_string = mod_options_string + mod_string;
		} else {
			mod_options_string = mod_options_string + ", " + mod_string;
		}
	}


	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute( "modification", xs_string, mod_options_string );

	attributes_for_parse_residue_selector( attlist, "residue_selector", "A residue selector specifying which residue or residues to predict on." );

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes(xsd, name_static(),
		"A metric for estimating the probability of a given site to be modified, as predicted by neural network.", attlist);
}
std::map< core::Size, core::Real >
PTMPredictionMetric::calculate( core::pose::Pose const & pose ) const {

	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );
	if ( residue_selector_ != nullptr ) {
		selection = residue_selector_->apply(pose);
	}
#ifdef USE_TENSORFLOW
	if ( modification_ == "Deamidation" ) {
		return tensorflow_protocol_->compute_deamidation_probability(pose, selection);
	} else {
		return tensorflow_protocol_->compute_ptm_probability(pose, selection, ptm_mode_ );
	}
#endif // USE_TENSORFLOW
	// define what happens if not compiled with tensorflow, need to make a dummy value to keep default compilation
	// howevere, I don't want the user to ever get that dummy map, so raise an exit instead
	utility_exit_with_message( "The PTMPredictionMetric needs to be compiled with either extras=tensorflow or extras=tensorflow_gpu!" );
	std::map< core::Size, core::Real > dummy_map;
	return dummy_map; // to keep compiler happy
}

/// @brief Set the residue selector that we'll be using, OP is used directly and the object is not cloned
/// @details Passing nullptr results in no residue selector being used.
void
PTMPredictionMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = selector_in;
}

/// @brief Get the residue selector.
/// @details If this returns nullptr, it means that no residue selector is being used.
core::select::residue_selector::ResidueSelectorCOP
PTMPredictionMetric::residue_selector() const {
	return residue_selector_;
}

/// @brief convert the type of modification enum to a string
std::string
PTMPredictionMetric::enum_to_string(
	PTMPredictionMetric::PTMPredictionMetricModes & ptm_mode
) {
	switch( ptm_mode ) {
	case PTMPredictionMetricModes::ACETYLATION :
		return "Acetylation";
	case PTMPredictionMetricModes::ARG_METHYLATION :
		return "ArgMethylation";
	case PTMPredictionMetricModes::CITRULLINATION :
		return "Citrullination";
	case PTMPredictionMetricModes::CROTONYLATION :
		return "Crotonylation";
	case PTMPredictionMetricModes::DEAMIDATION :
		return "Deamidation";
	case PTMPredictionMetricModes::GAMMA_CARBOXY_GLUTAMIC_ACID :
		return "GammaCarboxyGlutamicAcid";
	case PTMPredictionMetricModes::GLUTARYLATION :
		return "Glutarylation";
	case PTMPredictionMetricModes::GLUTATHIONYLATION :
		return "Glutathionylation";
	case PTMPredictionMetricModes::HYDROXYLATION :
		return "Hydroxylation";
	case PTMPredictionMetricModes::LYS_METHYLATION :
		return "LysMethylation";
	case PTMPredictionMetricModes::MALONYLATION :
		return "Malonylation";
	case PTMPredictionMetricModes::N_LINKED_GLYCOSYLATION :
		return "NlinkedGlycosylation";
	case PTMPredictionMetricModes::O_LINKED_GLYCOSYLATION :
		return "OlinkedGlycosylation";
	case PTMPredictionMetricModes::PHOSPHORYLATION :
		return "Phosphorylation";
	case PTMPredictionMetricModes::S_NITROSYLATION :
		return "SNitroysylation";
	case PTMPredictionMetricModes::SUCCINYLATION :
		return "Succinylation";
	case PTMPredictionMetricModes::SUMOYLATION :
		return "Sumoylation";
	case PTMPredictionMetricModes::UBIQUITINATION :
		return "Ubiquitination";
	case PTMPredictionMetricModes::INVALID_MODE :
		return "InvalidMode";
	}
	return ""; // to make compiler happy
}

PTMPredictionMetric::PTMPredictionMetricModes
PTMPredictionMetric::string_to_enum( std::string & modification_string ) {
	// loop through all possible prediction modes and check if equal
	for ( int i = 0; i <= static_cast< int >( PTMPredictionMetricModes::N_MODES ); i++ ) {
		PTMPredictionMetricModes mod = static_cast< PTMPredictionMetricModes >( i );
		if ( modification_string == enum_to_string( mod ) ) {
			return mod;
		}
	}

	utility_exit_with_message( "No valid modification type was selected!" );
	return PTMPredictionMetricModes::INVALID_MODE; // to keep compiler happy, should never be reached
}


void
PTMPredictionMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PTMPredictionMetric::provide_xml_schema( xsd );
}

std::string
PTMPredictionMetricCreator::keyname() const {
	return PTMPredictionMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PTMPredictionMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< PTMPredictionMetric >();
}

void
PTMPredictionMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"PTMPredictionMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the PTMPredictionMetric."
		)
	);
}
} //ptm_prediction
} //protocols

#ifdef    SERIALIZATION

template< class Archive >
void
protocols::ptm_prediction::PTMPredictionMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP ( modification_ ) );
	arc( CEREAL_NVP ( ptm_mode_ ) );
	arc( CEREAL_NVP (residue_selector_) );
	arc( CEREAL_NVP (tensorflow_protocol_) );

}

template< class Archive >
void
protocols::ptm_prediction::PTMPredictionMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( modification_ );
	arc( ptm_mode_ );
	core::select::residue_selector::ResidueSelectorOP residue_selector_local_;
	arc( residue_selector_local_);
	residue_selector_ = residue_selector_local_;
	arc( tensorflow_protocol_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::ptm_prediction::PTMPredictionMetric );
CEREAL_REGISTER_TYPE( protocols::ptm_prediction::PTMPredictionMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_ptm_prediction_PTMPredictionMetric )
#endif // SERIALIZATION





