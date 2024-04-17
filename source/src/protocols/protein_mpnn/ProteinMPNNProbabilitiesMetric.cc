// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/ProteinMPNNProbabilitiesMetric.cc
/// @brief A metric for predicting amino acid probabilities using ProteinMPNN (Dauparas et al.).
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// Unit headers
#include <protocols/protein_mpnn/ProteinMPNNProbabilitiesMetric.hh>
#include <protocols/protein_mpnn/ProteinMPNNProbabilitiesMetricCreator.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>
// Package Headers
#include <protocols/protein_mpnn/ProteinMPNN.hh>
#include <protocols/protein_mpnn/util.hh>
#include <protocols/protein_mpnn/ProteinMPNNMover.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/simple_moves/SimpleThreadingMover.hh>

// Utility Headers
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/Tracer.hh>
#include <core/pose/symmetry/util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility>
#include <utility/pointer/memory.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "protocols.protein_mpnn.ProteinMPNNProbabilitiesMetric" );


namespace protocols {
namespace protein_mpnn {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ProteinMPNNProbabilitiesMetric::ProteinMPNNProbabilitiesMetric():
	core::simple_metrics::PerResidueProbabilitiesMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ProteinMPNNProbabilitiesMetric::~ProteinMPNNProbabilitiesMetric()= default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ProteinMPNNProbabilitiesMetric::ProteinMPNNProbabilitiesMetric( ProteinMPNNProbabilitiesMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
ProteinMPNNProbabilitiesMetric::clone() const {
	return utility::pointer::make_shared< ProteinMPNNProbabilitiesMetric >( *this );
}

std::string
ProteinMPNNProbabilitiesMetric::name() const {
	return name_static();
}

std::string
ProteinMPNNProbabilitiesMetric::name_static() {
	return "ProteinMPNNProbabilitiesMetric";

}
std::string
ProteinMPNNProbabilitiesMetric::metric() const {
	return name_static();
}

void
ProteinMPNNProbabilitiesMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	SimpleMetric::parse_base_tag( tag );

	write_pssm_ = tag->getOption< std::string >( "write_pssm", "" );

	std::string residue_selector_rs_name = tag->getOption< std::string >( "residue_selector", "" );
	if ( residue_selector_rs_name.empty() ) {
		residue_selector_ = nullptr;
	} else {
		set_residue_selector( core::select::residue_selector::get_residue_selector( residue_selector_rs_name, datamap ));
	}

	std::string coord_selector_rs_name = tag->getOption< std::string >( "coord_selector", "" );
	if ( coord_selector_rs_name.empty() ) {
		coord_selector_rs_ = nullptr;
	} else {
		coord_selector_rs_ = core::select::residue_selector::get_residue_selector( coord_selector_rs_name, datamap );
	}

	std::string sequence_mask_selector_rs_name = tag->getOption< std::string >( "sequence_mask_selector", "" );
	if ( sequence_mask_selector_rs_name.empty() ) {
		sequence_mask_selector_rs_ = nullptr;
	} else {
		sequence_mask_selector_rs_ = core::select::residue_selector::get_residue_selector( sequence_mask_selector_rs_name, datamap );
	}

	tied_pos_rs_.clear();
	for ( utility::tag::TagCOP tied_def : tag->getTags( "TiedPositions" ) ) {
		if ( tied_def->hasOption( "residue_selectors" ) ) {
			utility::vector1< std::string > const res_selector_names =
				utility::string_split( tied_def->getOption< std::string >( "residue_selectors" ), ',' );

			utility::vector1< core::select::residue_selector::ResidueSelectorCOP > res_seles;
			for ( auto rs : res_selector_names ) {
				res_seles.emplace_back( core::select::residue_selector::get_residue_selector( rs, datamap ) );
			}

			tied_pos_rs_.emplace_back( res_seles );
		}
	}

}

void
ProteinMPNNProbabilitiesMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "write_pssm", xs_string, "Output filename for the psi-blast like position-specific-scoring-matrix to be used with the FavorSequenceProfile Mover", "" );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"residue_selector", xs_string,
		"A residue selector specifying which residue or residues to predict on. If none selected, all are passed.", "");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"coord_selector", xs_string,
		"Name of a residue selector that selects per-residue coordinates to pass to ProteinMPNN. If none selected, all coordinates are passed.", "");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"sequence_mask_selector", xs_string,
		"Name of a residue selector that selects positions to be masked.", "");

	// Attributes for TiedPositions
	AttributeList tied_pos_attrs;
	tied_pos_attrs + XMLSchemaAttribute( "residue_selectors", xs_string,
		"Comma separated list of residue selectors to tie together. "
		"The first residues of each selector will be tied together, then the second, etc. "
		"Each residue selector must have the same number of residues.");

	XMLSchemaSimpleSubelementList attrs_subelements;
	attrs_subelements.add_simple_subelement( "TiedPositions", tied_pos_attrs,
		"Used to define tied sets of residues." );
	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes_and_repeatable_subelements(xsd, name_static(),
		"A metric for estimating the probability of an amino acid at a given position, as predicted by ProteinMPNN.", attlist, attrs_subelements );

	//    XMLSchemaComplexTypeGenerator ct_gen;

	//    ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
	//            .element_name( name_static() )
	//            .description( "Author: Moritz Ertelt (moritz.ertelt@gmail.com)\n"
	//                          "A metric for estimating the probability of an amino acid at a given position, as predicted by the ProteinMPNN model." )
	//            .add_attributes( attlist )
	//            .add_optional_name_attribute()
	//            .set_subelements_repeatable( attrs_subelements, 0 )
	//            .write_complex_type_to_schema( xsd );
}

std::map< core::Size, std::map< core::chemical::AA, core::Real >>
ProteinMPNNProbabilitiesMetric::calculate( core::pose::Pose const &pose ) const {

	std::map< core::Size, std::map< core::chemical::AA, core::Real >> return_probabilities_map;

#ifdef USE_TORCH

    using namespace protocols::simple_moves;
    using namespace core::select;

    // define mpnn options
    ProteinMPNNOptions pose_options( pose );

    // set temperature to 1 so the probabilities don't get influenced
    pose_options.temperature = 1.0;

    // if specified by user, mask selected portions of the input sequence
    if ( sequence_mask_selector_rs_ != nullptr ) {
        utility::vector1< bool > sequence_mask = sequence_mask_selector_rs_->apply( pose );
        std::string masked_sequence;
        for ( core::Size position = 1; position <= pose.size(); ++position ) {
            if ( sequence_mask[position] ) {
                masked_sequence += 'X';
            } else {
                masked_sequence += pose.sequence()[position-1];
            }
        }
        pose_options.sequence = masked_sequence;
    }

    // use only specific backbone coordinates if specified by user
    if ( coord_selector_rs_ != nullptr ) {
        utility::vector1< bool > coord_mask = coord_selector_rs_->apply( pose );
        pose_options.coord_mask = coord_mask;
    }

    // convert tied positions as specified in tiedpositions tags to proteinmpnn tied positions
    utility::vector1< utility::vector1< core::Size > > tied_positions;

    for ( auto tied_res_selectors : tied_pos_rs_ ) {
        // getting indices from res selectors can occasionally be slow.
        // we store them here so we need to only compute them once
        utility::vector1 < utility::vector1< core::Size > > selection_positions = ProteinMPNNMover::residue_selectors_to_indices(pose, tied_res_selectors);

        core::Size rs_size = selection_positions[1].size();
        ProteinMPNNMover::validate_residue_selectors( selection_positions, rs_size );

        // for every residue of the selectors
        for ( core::Size residue_idx = 1; residue_idx <= rs_size; ++residue_idx ) {
            utility::vector1< core::Size > tied_resnums;

            // look at all the residue selectors in the tied_res_selectors and tie them together in order
            for ( auto pos_list : selection_positions ) {
                // put the residue number at index residue_idx of the residue selector into the tied_resnums
                tied_resnums.emplace_back( pos_list[ residue_idx ] );
            }

            tied_positions.emplace_back( tied_resnums );
        }
    }
    if ( core::pose::symmetry::is_symmetric( pose ) ) {
        core::conformation::symmetry::SymmetryInfoOP sym_info = core::pose::symmetry::symmetry_info( pose )->clone();

        // iterate over residues in asymmetric unit
        for ( core::Size ii = 1; ii <= sym_info->num_independent_residues(); ++ii ) {
            utility::vector1< core::Size > tied_set;
            tied_set.push_back( ii );

            // iterate over symmetric copies
            for ( core::Size subunit = 2; subunit <= sym_info->subunits(); ++subunit ) {
                core::Size clone = sym_info->equivalent_residue_on_subunit( subunit, ii );

                // tie symmetric positions
                tied_set.push_back( clone );
            }
            tied_positions.push_back( tied_set );
        }
    }

    pose_options.tied_positions = tied_positions;

    // predict the logits and probabilities (these will be for the complete pose)
    std::map<core::Size, utility::vector1<core::Real>> prob_map;
    std::map<core::Size, utility::vector1<core::Real>> logit_map;

    // ProteinMPNN inference
    ProteinMPNN::get_instance()->get_probabilities_and_logits( pose, pose_options, prob_map, logit_map );

    // remove positions if a selection is specified
    if ( residue_selector_ != nullptr ) {
        core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );
        selection = residue_selector_->apply(pose);
        for (auto it = prob_map.cbegin(); it != prob_map.cend();) {
            if (!selection[it->first]) {
                logit_map.erase(it->first);
                it = prob_map.erase(it);
            } else {
                ++it;
            }
        }
    }

    // write a PSSM in psi-blast format
    if ( !write_pssm_.empty() ) {
        // first get the sequence of all positions that were predicted on
        std::string selection_sequence = get_selection_sequence( pose.sequence(), logit_map );
        std::map< core::Size, std::map< core::chemical::AA, core::Real >> logit_aa_map;
        fill_return_map( logit_map, logit_aa_map );
        output_sequence_profile( selection_sequence, logit_aa_map, write_pssm_ );
    }

    // match the AA enums to the probabilities
    fill_return_map( prob_map, return_probabilities_map );

#else
	(void)pose; // avert error: unused parameter
	utility_exit_with_message( "This metric needs to be compiled with extras=pytorch!" );
#endif // USE_TORCH

	return return_probabilities_map;
}

#ifdef USE_TORCH
/// @brief Fill the return_map with the probabilities from the softmax_map for all amino acids
/// @param[in] softmax_map The map containing all probabilities or logits returned
/// @param[in] return_map The map that will be returned by the calculate function, will be filled with probabilities by this function
void
ProteinMPNNProbabilitiesMetric::fill_return_map( std::map< core::Size, utility::vector1< core::Real >> const &softmax_map, std::map< core::Size, std::map< core::chemical::AA, core::Real>> &return_map ) {
    // get the probabilities for all amino acids
    for ( auto const &pair: softmax_map) {
        core::Size selected_position = pair.first;
        utility::vector1<core::Real> const & softmax_vec = pair.second;
        for (core::Size prob_index = 0;
             prob_index < AA_ALPHABET.size() - 1; ++prob_index) { // -1 to avoid the unknown X type
            char const curr_aa = AA_ALPHABET[prob_index];
            core::chemical::AA const aa_enum = core::chemical::aa_from_oneletter_code(curr_aa);
            return_map[selected_position][aa_enum] = softmax_vec[prob_index + 1];
        }
    }
}
#endif // USE_TORCH

/// @brief Get the sequence of the selection from a map of positions->logits value
std::string
ProteinMPNNProbabilitiesMetric::get_selection_sequence( std::string const & pose_sequence, std::map< core::Size, utility::vector1< core::Real >> const & position_logit_map) {
	std::string selection_sequence;
	for ( auto const & pair : position_logit_map ) {
		selection_sequence += pose_sequence[pair.first - 1]; // pair.first is pose numbered position
	}
	return selection_sequence;
}

/// @brief Set the residue selector that we'll be using.
/// @details Passing nullptr results in no residue selector being used.
void
ProteinMPNNProbabilitiesMetric::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = std::move(selector_in);
}

void
ProteinMPNNProbabilitiesMetric::set_coord_selector_rs( core::select::residue_selector::ResidueSelectorCOP coord_selector_rs ) {
	coord_selector_rs_ = coord_selector_rs;
}

void
ProteinMPNNProbabilitiesMetric::set_sequence_mask_selector_rs( core::select::residue_selector::ResidueSelectorCOP sequence_mask_selector_rs ) {
	sequence_mask_selector_rs_ = sequence_mask_selector_rs;
}

void
ProteinMPNNProbabilitiesMetric::set_tied_pos_rs( utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs ) {
	tied_pos_rs_ = tied_pos_rs;
}

/// @brief Get the residue selector.
/// @details If this returns nullptr, it means that no residue selector is being used.
core::select::residue_selector::ResidueSelectorCOP
ProteinMPNNProbabilitiesMetric::residue_selector() const {
	return residue_selector_;
}

void
ProteinMPNNProbabilitiesMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ProteinMPNNProbabilitiesMetric::provide_xml_schema( xsd );
}

std::string
ProteinMPNNProbabilitiesMetricCreator::keyname() const {
	return ProteinMPNNProbabilitiesMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ProteinMPNNProbabilitiesMetricCreator::create_simple_metric() const {
	return utility::pointer::make_shared< ProteinMPNNProbabilitiesMetric >();
}

void
ProteinMPNNProbabilitiesMetric::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"ProteinMPNNProbabilitiesMetric", basic::citation_manager::CitedModuleType::SimpleMetric,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		""
		)
	);
#ifdef USE_TORCH
    citations.add( ProteinMPNN::get_ProteinMPNN_neural_net_citation() );
#endif //USE_TORCH
}
} //protein_mpnn
} //protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::protein_mpnn::ProteinMPNNProbabilitiesMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric>( this ) );

	arc( CEREAL_NVP ( write_pssm_ ) );
	arc( CEREAL_NVP ( residue_selector_ ) );
	arc( CEREAL_NVP ( sequence_mask_selector_rs_ ) );
	arc( CEREAL_NVP ( coord_selector_rs_ ) );
	arc( CEREAL_NVP ( tied_pos_rs_ ) );
}

template< class Archive >
void
protocols::protein_mpnn::ProteinMPNNProbabilitiesMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueProbabilitiesMetric >( this ) );

	arc( write_pssm_ );

	core::select::residue_selector::ResidueSelectorOP residue_selector_local_;
	arc( residue_selector_local_);
	residue_selector_ = residue_selector_local_;

	core::select::residue_selector::ResidueSelectorOP sequence_mask_selector_rs_local_;
	arc( sequence_mask_selector_rs_local_);
	sequence_mask_selector_rs_ = sequence_mask_selector_rs_local_;

	core::select::residue_selector::ResidueSelectorOP coord_selector_rs_local_;
	arc( coord_selector_rs_local_);
	coord_selector_rs_ = coord_selector_rs_local_;

	utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorOP > > tied_pos_rs_local_;
	arc( tied_pos_rs_local_ );
	tied_pos_rs_ = tied_pos_rs_local_;

}

SAVE_AND_LOAD_SERIALIZABLE( protocols::protein_mpnn::ProteinMPNNProbabilitiesMetric );
CEREAL_REGISTER_TYPE( protocols::protein_mpnn::ProteinMPNNProbabilitiesMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_protein_mpnn_ProteinMPNNProbabilitiesMetric )
#endif // SERIALIZATION





