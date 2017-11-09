// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/BondedResidueSelector.cc
/// @brief  Selects residues with chemical bonds to the passed input residues
/// @author Sharon Guffy (guffy@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/BondedResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/util.hh>
// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// C++ headers
#include <utility/assert.hh>


namespace core {
namespace select {
namespace residue_selector {

BondedResidueSelector::BondedResidueSelector():
	input_set_defined_(false)
{
}

BondedResidueSelector::BondedResidueSelector( std::set<core::Size> const & input_set)
{
	set_input_set(input_set);
}

BondedResidueSelector::BondedResidueSelector ( BondedResidueSelector const & other )
{
	if ( other.input_set_selector_defined() ) {
		use_input_set_selector_ = true;
		input_set_selector_ = other.input_set_selector();
	}
	input_set_defined_ = other.input_set_defined();
	set_input_set( other.input_set());
	input_set_str_ = other.input_set_string();

}


ResidueSelectorOP BondedResidueSelector::clone() const {
	return ResidueSelectorOP ( new BondedResidueSelector( *this ) );
}


BondedResidueSelector::~BondedResidueSelector() {}

ResidueSubset
BondedResidueSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( input_set_defined_ );

	ResidueSubset subset( pose.total_residue(), false );
	//This is in a separate variable in case we're getting input_set from a ResidueSelector
	std::set< Size > input_set_tmp;

	// set subset to input_set.
	get_input_set(pose, subset, input_set_tmp);
	//Get all residues bonded to residues in input_set_
	for ( std::set< core::Size>::const_iterator it = input_set_tmp.begin(); it != input_set_tmp.end(); ++it ) {
		//First we need to get all the AtomIDs for this residue
		//Format = (atomno, resno)
		for ( core::Size atomnum = 1; atomnum <= pose.residue(*it).natoms(); ++atomnum ) {

			//This method returns an AtomID
			utility::vector1< core::id::AtomID > bonded_atoms = pose.conformation().bonded_neighbor_all_res( core::id::AtomID( atomnum, *it ) );
			//This will probably end up being redundant, but we're just changing a boolean value, so that's okay
			for ( core::Size i = 1; i <= bonded_atoms.size(); ++i ) {
				subset[ bonded_atoms[ i ].rsd() ] = true;
			}
		}
	}


	return subset;
}

void
BondedResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( tag->hasOption("residue_selector") ) {
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "BondedResidueSelector takes EITHER 'selector' OR 'resnum' options, not both!\n" );
		}
		if ( tag->size() > 1 ) { // 1 if no subtags exist
			throw utility::excn::EXCN_Msg_Exception( "BondedResidueSelector can only have one ResidueSelector loaded!\n" );
		}
		// grab the ResidueSelector from the selector option
		// and then grab each of the indicated residue selectors from the datamap.
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "residue_selector" );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selector' from BondedResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}

		try {
			ResidueSelectorCOP selector = datamap.get_ptr< ResidueSelector const >( "ResidueSelector", selector_str );
			set_input_set_selector( selector );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from BondedResidueSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // get focus selector from tag
		if ( tag->hasOption("resnums") ) {
			throw utility::excn::EXCN_Msg_Exception( "BondedResidueSelector takes EITHER a 'resnums' tag or a selector subtag, not both!\n" );
		}

		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw utility::excn::EXCN_Msg_Exception( "BondedResidueSelector takes at most one ResidueSelector to determine the input_set!\n" );
		}
		ResidueSelectorCOP rs = ResidueSelectorFactory::get_instance()->new_residue_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_input_set_selector( rs );

	} else { // do not get input_set from ResidueSelectors but load resnums string instead
		try {
			set_input_set ( tag->getOption< std::string >( "resnums" ) );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream err_msg;
			err_msg << "Failed to access option 'resnums' from BondedResidueSelector::parse_my_tag.\n";
			err_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
	}

}
///@brief Sets subset to be true for values contained in input_set
void
BondedResidueSelector::get_input_set(
	core::pose::Pose const & pose,
	ResidueSubset & subset,
	std::set< Size > & input_set
) const
{
	if ( input_set_selector_ && use_input_set_selector_ ) {
		subset = input_set_selector_->apply( pose );

		for ( Size ii = 1; ii <= subset.size(); ++ii ) {
			if ( subset[ ii ] ) {
				input_set.insert( ii );
			}
		}
	} else { // grab from set and string
		std::set< Size > const res_vec( get_resnum_list( input_set_str_, pose ) );
		input_set.insert( input_set_.begin(), input_set_.end() );
		input_set.insert( res_vec.begin(), res_vec.end() );
		for ( std::set< Size >::const_iterator it = input_set.begin();
				it != input_set.end(); ++it ) {
			if ( *it == 0 || *it > subset.size() ) {
				std::stringstream err_msg;
				err_msg << "Residue " << *it << " not found in pose!\n";
				throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
			}
			subset[ *it ] = true; // may want to use a tmp subset so we don't wind up with a half-set subset
		}
	}
}

void
BondedResidueSelector::set_input_set( std::set<Size> const & input_set )
{
	input_set_ = input_set;
	input_set_defined_ = true;
	use_input_set_selector_ = false;
}

void
BondedResidueSelector::set_input_set( std::string const & input_set_str )
{
	input_set_str_ = input_set_str;
	input_set_defined_ = true;
	use_input_set_selector_ = false;
}

void BondedResidueSelector::set_input_set_selector( ResidueSelectorCOP rs )
{
	input_set_selector_ = rs;
	input_set_defined_ = true;
	use_input_set_selector_ = true;
}


std::string BondedResidueSelector::get_name() const {
	return BondedResidueSelector::class_name();
}

std::string BondedResidueSelector::class_name() {
	return "Bonded";
}



ResidueSelectorOP
BondedResidueSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new BondedResidueSelector );
}

std::string
BondedResidueSelectorCreator::keyname() const {
	return BondedResidueSelector::class_name();
}

std::set< core::Size > BondedResidueSelector::input_set() const
{
	return input_set_;
}
core::select::residue_selector::ResidueSelectorCOP BondedResidueSelector::input_set_selector() const
{
	return input_set_selector_;
}
std::string BondedResidueSelector::input_set_string() const
{
	return input_set_str_;
}
bool BondedResidueSelector::input_set_defined() const
{
	return input_set_defined_;
}
bool BondedResidueSelector::input_set_selector_defined() const
{
	return use_input_set_selector_;
}

void
BondedResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "residue_selector", xs_string, "Name of residue selector specifying residues for which to select bonded partners" );
	core::pose::attributes_for_get_resnum_selector( attributes, xsd, "resnums" );
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(), "Selects all residues that are attached to the provided residues by a chemical bond", attributes );

}



void
BondedResidueSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	BondedResidueSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core
