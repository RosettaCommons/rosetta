// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/residue_selectors/HBondSelector.hh
/// @brief  A ResidueSelector that selects residues forming hydrogen bonds with specified residues. By default it ignores backbone-backbone hydrogen bonds.
/// @author Sharon Guffy (guffy@email.unc.edu)



// Unit headers
// Unit headers
#include <protocols/residue_selectors/HBondSelector.hh>
#include <protocols/residue_selectors/HBondSelectorCreator.hh>


//Protocols
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <core/types.hh>
// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace residue_selectors {


HBondSelector::HBondSelector():
	hbond_energy_cutoff_( -0.5 ),
	include_bb_bb_( false ),
	input_set_defined_( false ),
	use_input_set_selector_( true )
{
	//Set default score function
	scorefxn_ = core::scoring::get_score_function();
	//If the input set is not defined, then we'll use all residues
	input_set_selector_ = core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector );
}

HBondSelector::HBondSelector( HBondSelector const & src ){
	//Set all data members from values in src
	set_hbond_energy_cutoff( src.get_hbond_energy_cutoff() );
	set_include_bb_bb( src.get_include_bb_bb() );
	set_scorefxn( src.get_scorefxn() );
	set_input_set_str( src.get_input_set_str() );
	set_input_set_selector( src.get_input_set_selector() );
	input_set_defined_ = src.get_input_set_defined();
	use_input_set_selector_ = src.get_use_input_set_selector();
}

HBondSelector::~HBondSelector(){}

core::select::residue_selector::ResidueSelectorOP
HBondSelector::clone() const{
	return core::select::residue_selector::ResidueSelectorOP( new HBondSelector( *this ) );
}

core::select::residue_selector::ResidueSubset
HBondSelector::apply( core::pose::Pose const & pose ) const{
	//Initialize input set
	core::pose::PoseOP copy_pose = pose.clone();
	std::set< Size > input_set;
	compute_input_set( pose, input_set );

	//Make sure the score function is correctly initialized
	//Sorry, you're on your own!
	//core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, scorefxn_ );
	core::scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
	energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
	scorefxn_->set_energy_method_options( energymethodoptions );
	scorefxn_->score( *copy_pose );
	//Now we'll find the hydrogen bonds
	core::scoring::hbonds::HBondSet hbond_set;
	hbond_set.setup_for_residue_pair_energies( *copy_pose, false, false );

	//Initialize return value
	core::select::residue_selector::ResidueSubset subset( pose.total_residue(), false );

	//Iterate over hbonds in the set
	for ( core::Size i = 1; i <= hbond_set.nhbonds(); ++i ) {
		core::scoring::hbonds::HBond const & hbond = hbond_set.hbond(i);
		//Ignore bb-bb if appropriate
		if ( !include_bb_bb_ && hbond.don_hatm_is_backbone() && hbond.acc_atm_is_backbone() ) {
			continue;
		}
		//Check that the energy is low enough
		if ( hbond.energy() > hbond_energy_cutoff_ ) {
			continue;
		}
		//Otherwise count anything involving residues in input_set
		if ( input_set.find( hbond.acc_res() ) != input_set.end() ) {
			subset[ hbond.don_res() ] = true;
		}
		if ( input_set.find( hbond.don_res() ) != input_set.end() ) {
			subset[ hbond.acc_res() ] = true;
		}
	}
	return subset;
}

void
HBondSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	if ( tag->hasOption("residue_selector") ) {
		if ( tag->hasOption("resnums") ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "HBondResidueSelector takes EITHER 'selector' OR 'resnum' options, not both!\n" );
		}
		if ( tag->size() > 1 ) { // 1 if no subtags exist
			throw CREATE_EXCEPTION(utility::excn::Exception,  "HBondResidueSelector can only have one ResidueSelector loaded!\n" );
		}
		// grab the ResidueSelector from the selector option
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "residue_selector" );
		} catch ( utility::excn::Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selector' from HBondResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
		try {
			core::select::residue_selector::ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
			set_input_set_selector( selector );
		} catch ( utility::excn::Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from HBondResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // get input selector from tag
		if ( tag->hasOption("resnums") ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "HBondResidueSelector takes EITHER a 'resnums' tag or a selector subtag, not both!\n" );
		}
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "HBondResidueSelector takes at most one ResidueSelector to determine the input_set!\n" );
		}
		core::select::residue_selector::ResidueSelectorCOP rs = core::select::residue_selector::ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_input_set_selector( rs );
	} else if ( tag->hasOption( "resnums" ) ) { // do not get input_set from ResidueSelectors but load resnums string instead
		try {
			set_input_set_str ( tag->getOption< std::string >( "resnums" ) );
		} catch ( utility::excn::Exception e ) {
			std::stringstream err_msg;
			err_msg << "Failed to access option 'resnums' from HBondResidueSelector::parse_my_tag.\n";
			err_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg.str() );
		}
	}
	//Set other options
	//
	hbond_energy_cutoff_ = tag->getOption< core::Real >( "hbond_energy_cutoff", -0.5 );
	include_bb_bb_ = tag->getOption< bool >( "include_bb_bb", false );

	//Score Function
	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );

	//If the input set has not been defined at this point, then we'll use all residues
	if ( !input_set_defined_ ) {
		//Okay, so I can't actually do this in the apply method--what now?
		use_input_set_selector_ = true;
		input_set_selector_ = core::select::residue_selector::ResidueSelectorCOP( new core::select::residue_selector::TrueResidueSelector );
	}

}
///@details Non-static function to return class name
std::string
HBondSelector::get_name() const{
	return HBondSelector::class_name();
}

/// @details Static function to return class name; anything else returning the class name should call this function to avoid discrepancies
std::string
HBondSelector::class_name(){
	return "HBond";
}

/// @details Static function allowing evaluation of XML before parsing
void
HBondSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "include_bb_bb", xsct_rosetta_bool, "Include backbone-backbone hydrogen bonds?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "hbond_energy_cutoff", xsct_real, "Upper energy cutoff for whether we will count a hydrogen bond", "-0.5" )
		+ XMLSchemaAttribute( "residue_selector", xs_string, "Name of residue selector specifying residues for which to select hydrogen bonded partners" );
	core::pose::attributes_for_get_resnum_selector( attlist, xsd, "resnums" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	core::select::residue_selector::xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name() ,"A ResidueSelector that selects all residues with hydrogen bonds to the specified starting residues. If no input residues are specified, selects all residues forming hydrogen bonds.", attlist );

}

//Getters and Setters
void
HBondSelector::set_hbond_energy_cutoff( core::Real input_value ){
	hbond_energy_cutoff_ = input_value;
}

core::Real
HBondSelector::get_hbond_energy_cutoff() const{
	return hbond_energy_cutoff_;
}

void
HBondSelector::set_include_bb_bb ( bool const input_setting ){
	include_bb_bb_ = input_setting;
}

bool
HBondSelector::get_include_bb_bb() const{
	return include_bb_bb_;
}

void
HBondSelector::set_scorefxn( core::scoring::ScoreFunctionCOP sfxn ){
	scorefxn_ = sfxn->clone();
}

core::scoring::ScoreFunctionCOP
HBondSelector::get_scorefxn() const{
	return scorefxn_;
}

std::string
HBondSelector::get_input_set_str() const{
	return input_set_str_;
}

void
HBondSelector::set_input_set_str( std::string input ){
	input_set_str_ = input;
	input_set_defined_ = true;
	use_input_set_selector_ = false;
}

core::select::residue_selector::ResidueSelectorCOP
HBondSelector::get_input_set_selector() const{
	return input_set_selector_;
}

void
HBondSelector::set_input_set_selector( core::select::residue_selector::ResidueSelectorCOP select ){
	input_set_selector_ = select;
	input_set_defined_ = true;
	use_input_set_selector_ = true;
}
//We need these two getters for the copy constructor only (they shouldn't be set externally)
bool
HBondSelector::get_input_set_defined() const{
	return input_set_defined_;
}

bool
HBondSelector::get_use_input_set_selector() const{
	return use_input_set_selector_;
}
void
HBondSelector::compute_input_set( core::pose::Pose const & pose, std::set< core::Size > & input_set ) const{
	if ( input_set_selector_ && use_input_set_selector_ ) {
		core::select::residue_selector::ResidueSubset subset = input_set_selector_->apply( pose );
		for ( core::Size ii = 1; ii <= subset.size(); ++ii ) {
			if ( subset[ ii ] ) {
				input_set.insert( ii );
			}
		}
	} else { // grab from string
		std::set< Size > const string_set( get_resnum_list( input_set_str_, pose ) );
		input_set.insert( string_set.begin(), string_set.end() );

		for ( core::Size res: input_set ) {
			if ( res == 0 || res > pose.total_residue() ) {
				std::stringstream err_msg;
				err_msg << "Residue " << res << " not found in pose!\n";
				throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg.str() );
			}
		}
	}
}
//CREATOR METHODS
core::select::residue_selector::ResidueSelectorOP
HBondSelectorCreator::create_residue_selector() const{
	return core::select::residue_selector::ResidueSelectorOP( new HBondSelector );
}

std::string
HBondSelectorCreator::keyname() const{
	return HBondSelector::class_name();
}

void
HBondSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	HBondSelector::provide_xml_schema( xsd );
}

} //protocols
} //residue_selectors


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::residue_selectors::HBondSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( hbond_energy_cutoff_ ) ); // core::Real
	arc( CEREAL_NVP( include_bb_bb_ ) ); // bool
	arc( CEREAL_NVP( scorefxn_ ) ); // core::scoring::ScoreFunctionOP
	arc( CEREAL_NVP( input_set_selector_ ) ); //core::select::residue_selector::ResidueSelectorCOP
	arc( CEREAL_NVP( input_set_str_ ) ); //std::string
	arc( CEREAL_NVP( input_set_defined_ ) ); // bool
	arc( CEREAL_NVP( use_input_set_selector_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::residue_selectors::HBondSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( hbond_energy_cutoff_ ); // core::Real
	arc( include_bb_bb_ ); // bool
	arc( scorefxn_ ); // core::scoring::ScoreFunctionOP
	arc( input_set_selector_ ); //core::select::residue_selector::ResidueSelectorCOP
	arc( input_set_str_ ); //std::string
	arc( input_set_defined_ ); // bool
	arc( use_input_set_selector_ ); // bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::residue_selectors::HBondSelector );
CEREAL_REGISTER_TYPE( protocols::residue_selectors::HBondSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_residue_selectors_HBondSelector )
#endif // SERIALIZATION
