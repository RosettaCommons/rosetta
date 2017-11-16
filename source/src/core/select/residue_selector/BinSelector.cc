// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/BinSelector.hh
/// @brief  A ResidueSelector that selects residues based on their torsion bin (e.g. ABEGO bin).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <core/select/residue_selector/BinSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.select.residue_selector.BinSelector" );


namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

/// @brief Constructor.
///
BinSelector::BinSelector() :
	initialized_(false),
	bin_transition_calculator_( new core::scoring::bin_transitions::BinTransitionCalculator() ),
	select_only_alpha_aas_(true),
	bin_name_(""),
	bin_params_file_name_("ABEGO")
	//TODO -- initialize all vars here.
{}

/// @brief Copy constructor.
///
BinSelector::BinSelector( BinSelector const &src ) :
	initialized_( src.initialized_ ),
	bin_transition_calculator_( src.bin_transition_calculator_->clone() ),
	select_only_alpha_aas_(src.select_only_alpha_aas_),
	bin_name_( src.bin_name_ ),
	bin_params_file_name_( src.bin_params_file_name_ )
	//TODO -- initialize all vars here.
{}

/// @brief Destructor.
///
BinSelector::~BinSelector() {}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
BinSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast<core::select::residue_selector::ResidueSelector>(
		BinSelectorOP( new BinSelector(*this) )
		)
	);
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ResidueSubset
BinSelector::apply(
	core::pose::Pose const & pose
) const {
	runtime_assert_string_msg( initialized() && bin_transition_calculator_->bin_params_loaded(),
		"Error in core::select::residue_selector::BinSelector::apply(): The parameters for the BinTransitionCalculator have not been loaded." );
	core::Size const nres( pose.size() ); //Number of residues in the pose
	ResidueSubset outvect( nres, false ); //Output vector, initialized to a vector of "false".
	for ( core::Size i=1; i<=nres; ++i ) {
		// Ignore ligands, non-alpha AAs (if instructed), residues without a
		// LOWER or UPPER, and residues with broken LOWER or UPPER.
		if ( !pose.residue_type( i ).is_polymer() ) continue;
		if ( select_only_alpha_aas() && !pose.residue_type(i).is_alpha_aa() ) continue;
		if ( !pose.residue(i).has_lower_connect() || !pose.residue(i).has_upper_connect() ) continue;
		if ( !pose.residue(i).connected_residue_at_resconn( pose.residue_type(i).lower_connect_id() ) ) continue;
		if ( !pose.residue(i).connected_residue_at_resconn( pose.residue_type(i).upper_connect_id() ) ) continue;

		core::Size data(0), bin(0);
		bin_transition_calculator_->find_data_and_bin( bin_name(), pose.residue(i), data, bin, false);
		if ( !data || !bin ) continue; //Do not select this residue if no bin exists for this type.

		//Set the selection based on whether this residue is in the bin:
		outvect[i] = bin_transition_calculator_->is_in_bin( pose.residue(i), data, bin, false );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "BinSelector has selected:" << std::endl;
		for ( core::Size i=1, imax=outvect.size(); i<=imax; ++i ) {
			TR.Debug << i << "\t" << (outvect[i] ? "TRUE" : "FALSE");
			if ( pose.residue_type(i).is_alpha_aa() ) {
				TR.Debug << "\tphi:" << pose.phi(i) << "\tpsi:" << pose.psi(i) << "\tomega:" << pose.omega(i) << std::endl;
			}
			TR.Debug << std::endl;
		}
		TR.Debug.flush();
	}

	return outvect;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
BinSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*datamap*/
) {
	//Get the bin params file:
	set_bin_params_file_name( tag->getOption< std::string >( "bin_params_file", bin_params_file_name() ) );

	//Set desired bin.
	set_bin_name( tag->getOption<std::string>("bin") );

	//Set whether we're only selecting alpha-amino acids.
	set_select_only_alpha_aas( tag->getOption< bool >( "select_only_alpha_aas", select_only_alpha_aas() ) );

	//Actually load the bin params file.
	initialize_and_check();
}

/// @brief Load the bin params file and check that settings are consistent.
/// @details Must be called before apply() function.
void
BinSelector::initialize_and_check() {
	if ( initialized() || bin_transition_calculator_->bin_params_loaded() ) {
		utility_exit_with_message(
			"Error in core::select::residue_selector::BinSelector::initialize_and_check(): The bin params file was already loaded!  This operation cannot be repeated.");
	}
	bin_transition_calculator_->load_bin_params( bin_params_file_name() );

	runtime_assert_string_msg( bin_transition_calculator_->bin_definition_exists( bin_name() ),
		"Error in core::select::residue_selector::BinSelector::initialize_and_check(): The bin params file defines no bin named " + bin_name() + "." );

	if ( TR.visible() ) {
		TR << "Loaded bin parameters file \"" << bin_params_file_name() << "\" and set bin to select to \"" << bin_name() << "\".";
		if ( select_only_alpha_aas() ) {
			TR << "  Only selecting alpha-amino acids." << std::endl;
		} else {
			TR << "  Selecting all polymeric types in the specified bin." << std::endl;
		}
	}

	//Set that we've initialized this BinSelector.
	initialized_=true;

	return;
}

/// @brief Load the bin params file baed on a file contents string (instead of loading directly
/// from disk) and check that settings are consistent.
/// @details Must be called as an alternative to initialize_and_check() before apply() function.
void
BinSelector::initialize_from_file_contents_and_check (
	std::string const &filecontents
) {
	if ( initialized() || bin_transition_calculator_->bin_params_loaded() ) {
		utility_exit_with_message(
			"Error in core::select::residue_selector::BinSelector::initialize_from_file_contents_and_check(): The bin params file was already loaded!  This operation cannot be repeated.");
	}
	bin_transition_calculator_->load_bin_params_from_file_contents( filecontents );

	runtime_assert_string_msg( bin_transition_calculator_->bin_definition_exists( bin_name() ),
		"Error in core::select::residue_selector::BinSelector::initialize_from_file_contents_and_check(): The bin params file defines no bin named " + bin_name() + "." );

	if ( TR.visible() ) {
		TR << "Loaded bin parameters from file contents and set bin to select to \"" << bin_name() << "\".";
		if ( select_only_alpha_aas() ) {
			TR << "  Only selecting alpha-amino acids." << std::endl;
		} else {
			TR << "  Selecting all polymeric types in the specified bin." << std::endl;
		}
	}

	//Set that we've initialized this BinSelector.
	initialized_=true;

	return;
}

/// @brief Get the mover class name.
///
std::string BinSelector::get_name() const {
	return BinSelector::class_name();
}

/// @brief Get the mover class name.
///
std::string BinSelector::class_name() {
	return "Bin";
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
BinSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::attribute_w_default(
		"bin_params_file", xs_string,
		"The filename of a bin_params file that defines a set of torsion bins. "
		"Current options include \"ABEGO\", \"ABBA\" (a symmetric version of "
		"the ABEGO nomenclature useful for mixed D/L design), and \"PRO_DPRO\" "
		"(which defines bins \"LPRO\" and \"DPRO\" corresponding to the regions "
		"of Ramachandran space accessible to L- and D-proline, respectively). "
		"Predefined bin_params files are in database/protocol_data/generalizedKIC/bin_params/. "
		"A custom-written bin_params file may also be provided with its relative "
		"path from the execution directory.",
		"ABEGO"  ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default(
		"bin", xs_string,
		"The name of the mainchain torsion bin.",
		""));
	attributes.push_back( XMLSchemaAttribute::attribute_w_default(
		"select_only_alpha_aas", xsct_rosetta_bool,
		"If true, only alpha-amino acids are selected. If false, "
		"any polymeric type allowed by the bin definitions file is selected.",
		"true"  ));
	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The BinSelector selects residues that fall in a named mainchain torsion bin "
		"(e.g. the \"A\" bin, corresponding to alpha-helical residues by the \"ABEGO\" "
		"nomenclature). Non-polymeric residues are ignored. By default, only "
		"alpha-amino acids are selected, though this can be disabled.",
		attributes );
}

core::select::residue_selector::ResidueSelectorOP
BinSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP(
		utility::pointer::dynamic_pointer_cast< core::select::residue_selector::ResidueSelector > (
		BinSelectorOP( new BinSelector )
		)
	);
}

std::string
BinSelectorCreator::keyname() const {
	return BinSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
///
void
BinSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	BinSelector::provide_xml_schema( xsd );
}


} //core
} //select
} //residue_selector

#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::BinSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( initialized_ ) ); // bool
	arc(  CEREAL_NVP( bin_transition_calculator_ ) ); // core::scoring::bin_transitions::BinTransitionCalculatorOP
	arc( CEREAL_NVP( select_only_alpha_aas_ ) ); //bool
	arc( CEREAL_NVP( bin_name_ ) ); //std::string
	arc( CEREAL_NVP( bin_params_file_name_ ) ); //std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::BinSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( initialized_ ); // bool
	arc( bin_transition_calculator_ ); // core::scoring::bin_transitions::BinTransitionCalculatorOP
	arc( select_only_alpha_aas_ ); //bool
	arc( bin_name_ ); //std::string
	arc( bin_params_file_name_ ); //std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::BinSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::BinSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_BinSelector )
#endif // SERIALIZATION
