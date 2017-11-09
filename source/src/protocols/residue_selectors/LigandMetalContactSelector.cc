// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/residue_selectors/LigandMetalContactSelector.hh
/// @brief  This residue selector takes a selector or residue number of a ligand and returns any residues in contact with metal atoms in the ligand.
/// @author Allison Watwood (acw538@msstate.edu)

// Unit headers
#include <protocols/residue_selectors/LigandMetalContactSelector.hh>
#include <protocols/residue_selectors/LigandMetalContactSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>
#include <core/util/metalloproteins_util.hh>
#include <core/pose/selection.hh>
#include <core/id/AtomID.hh>

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


static THREAD_LOCAL basic::Tracer TR( "protocols.residue_selectors.LigandMetalContactSelector" );

namespace protocols {
namespace residue_selectors {

/// @brief Constructor.
///
LigandMetalContactSelector::LigandMetalContactSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
LigandMetalContactSelector::~LigandMetalContactSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
LigandMetalContactSelector::LigandMetalContactSelector(LigandMetalContactSelector const & src):
	core::select::residue_selector::ResidueSelector( src )
{
	using_a_resselect_ = src.using_a_resselect_;
	dist_cutoff_multiplier_  = src.dist_cutoff_multiplier_;
	input_set_selector_ = src.input_set_selector_;
	resnum_string_ = src.resnum_string_;

}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
LigandMetalContactSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new LigandMetalContactSelector(*this) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
LigandMetalContactSelector::ResidueSubset
LigandMetalContactSelector::apply(core::pose::Pose const & pose) const {

	ResidueSubset subset(pose.total_residue(), false);
	std::set<Size> input_set_tmp;
	calculate_ligand_resnums(pose, input_set_tmp);

	for ( std::set<core::Size>::const_iterator it = input_set_tmp.begin(); it != input_set_tmp.end(); ++it ) {

		//for(core::Size atomnum = 1; atomnum <= pose.residue(*it).natoms(); ++atomnum){
		std::map<core::Size, utility::vector1<core::id::AtomID>> contact_map = core::util::find_metalbinding_atoms_for_complex(pose, *it, dist_cutoff_multiplier_);

		for ( std::pair<core::Size, utility::vector1<core::id::AtomID>> contact_pair:contact_map ) {

			for ( core::Size i = 1; i <= contact_pair.second.size(); ++i ) {

				subset[contact_pair.second[i].rsd()] = true;
				TR << "Residue " << i << " selected"<< std::endl;
			}
		}
	}
	return subset;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
LigandMetalContactSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	if ( tag->hasOption("residue_selector") ) {
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "LigandMetalContactSelector takes EITHER 'selector' OR 'resnums' options, not both!\n" );
		}
		if ( tag->size() > 1 ) { // 1 if no subtags exist
			throw utility::excn::EXCN_Msg_Exception( "LigandMetalContactSelector can only have one ResidueSelector loaded!\n" );
		}
		// grab the ResidueSelector from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str = "";
		try {
			selector_str = tag->getOption< std::string >( "residue_selector" );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selector' from LigandMetalContactSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}

		try {
			core::select::residue_selector::ResidueSelectorCOP selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_str );
			set_input_set_selector( selector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from LigandMetalContactSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // get focus selector from tag
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "LigandMetalContactSelector takes EITHER a 'resnums' tag or a selector subtag, not both!\n" );
		}

		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw utility::excn::EXCN_Msg_Exception( "LigandMetalContactSelector takes at most one ResidueSelector to determine the input_set!\n" );
		}
		core::select::residue_selector::ResidueSelectorCOP rs = core::select::residue_selector::ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_input_set_selector( rs );
	} else { // do not get input_set from ResidueSelectors but load resnums string instead
		try {
			set_resnum_string ( tag->getOption< std::string >( "resnums" ) );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream err_msg;
			err_msg << "Failed to access option 'resnums' from LigandMetalContactSelector::parse_my_tag.\n";
			err_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
	}

	dist_cutoff_multiplier_ = tag->getOption< core::Real >( "dist_cutoff_multiplier", 1.0 );

}

std::string LigandMetalContactSelector::get_name() const
{
	return LigandMetalContactSelector::class_name();
}

std::string LigandMetalContactSelector::class_name()
{
	return "LigandMetalContactSelector";
}

void LigandMetalContactSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "residue_selector", xs_string, "Name of the residue selector for the ligand")
		+ XMLSchemaAttribute::attribute_w_default( "dist_cutoff_multiplier", xsct_real, "Multiplier for the distance from the metal atom for contact detection", "1");
	core::pose::attributes_for_get_resnum_selector( attributes, xsd, "resnums");
	core::select::residue_selector::xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), "This residue selector selects for the residues in contact with the ligand metal", attributes );


}

core::select::residue_selector::ResidueSelectorOP
LigandMetalContactSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new LigandMetalContactSelector );
}

std::string
LigandMetalContactSelectorCreator::keyname() const {
	return LigandMetalContactSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
LigandMetalContactSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	LigandMetalContactSelector::provide_xml_schema( xsd );
}

void
LigandMetalContactSelector::calculate_ligand_resnums(core::pose::Pose const & pose, std::set<Size> & input_set) const{
	ResidueSubset local_subset(pose.total_residue(), false);

	if ( input_set_selector_ && using_a_resselect_ ) {
		local_subset = input_set_selector_->apply(pose);

		for ( Size ii = 1; ii <= local_subset.size(); ++ii ) {
			if ( local_subset[ii] ) {
				input_set.insert(ii);
			}
		}
	} else {
		std::set<Size> const res_vec( core::pose::get_resnum_list( resnum_string_, pose));
		input_set.insert(res_vec.begin(), res_vec.end());
		for ( std::set<Size>::const_iterator it = input_set.begin(); it != input_set.end(); ++it ) {
			if ( *it == 0 || *it > local_subset.size() ) {
				std::stringstream err_msg;
				err_msg << "Residue " << *it <<" not found in pose!/n";
				throw utility::excn::EXCN_Msg_Exception(err_msg.str());
			}
		}
	}
}


std::string
LigandMetalContactSelector::get_resnum_string() const{
	return resnum_string_;
}

void
LigandMetalContactSelector::set_resnum_string(std::string const & input_string){
	resnum_string_ = input_string;
	using_a_resselect_ = false;
}

core::select::residue_selector::ResidueSelectorCOP
LigandMetalContactSelector::get_input_set_selector() const{
	return input_set_selector_;
}

void
LigandMetalContactSelector::set_input_set_selector( core::select::residue_selector::ResidueSelectorCOP rs){
	input_set_selector_ = rs;
	using_a_resselect_ = true;
}

bool
LigandMetalContactSelector::get_using_a_resselect() const{
	return using_a_resselect_;
}

void
LigandMetalContactSelector::set_using_a_resselect(bool new_bool){
	using_a_resselect_ = new_bool;
}

core::Real
LigandMetalContactSelector::get_dist_cutoff_multiplier() const{
	return dist_cutoff_multiplier_;
}

void
LigandMetalContactSelector::set_dist_cutoff_multiplier(core::Real const new_d_c_m){
	dist_cutoff_multiplier_ = new_d_c_m;
}

} //protocols
} //residue_selectors

#ifdef    SERIALIZATION

// See the serialization documentation here.  There is a script you can run.
//  https://wiki.rosettacommons.org/index.php/SerializationFAQ

template< class Archive >
void
protocols::residue_selectors::LigandMetalContactSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	// You need to add "arc( CEREAL_NVP( datamember_name ) );" calls here for each of your data members.
	arc( CEREAL_NVP( resnum_string_ )); //std::string
	arc( CEREAL_NVP( input_set_selector_ )); //ResidueSelectorCOP
	arc( CEREAL_NVP( using_a_resselect_ )); //bool
	arc( CEREAL_NVP( dist_cutoff_multiplier_ )); //core::Real
}

template< class Archive >
void
protocols::residue_selectors::LigandMetalContactSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	// You need to add "arc( datamember_name );" calls here for each of your data members.
	arc( resnum_string_ ); //std::string
	arc( input_set_selector_ ); //ResidueSelectorCOP
	arc( using_a_resselect_ ); //bool
	arc( dist_cutoff_multiplier_ ); //core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::residue_selectors::LigandMetalContactSelector );
CEREAL_REGISTER_TYPE( protocols::residue_selectors::LigandMetalContactSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_residue_selectors_LigandMetalContactSelector )
#endif // SERIALIZATION
