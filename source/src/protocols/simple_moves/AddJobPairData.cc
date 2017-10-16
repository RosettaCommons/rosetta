// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddJobPairData.cc
/// @brief  A really simple mover that takes some data in through xml and appends it to the pose
/// @author Sam DeLuca <Samuel.l.deluca@vanderbilt.edu)

#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <protocols/simple_moves/AddJobPairData.hh>
#include <protocols/simple_moves/AddJobPairDataCreator.hh>
#include <protocols/jd2/util.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.AddJobPairData" );

// XRW TEMP std::string AddJobPairDataCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AddJobPairData::mover_name();
// XRW TEMP }

// XRW TEMP moves::MoverOP AddJobPairDataCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return moves::MoverOP( new AddJobPairData );
// XRW TEMP }

// XRW TEMP std::string AddJobPairData::mover_name()
// XRW TEMP {
// XRW TEMP  return "AddJobPairData";
// XRW TEMP }

AddJobPairData::AddJobPairData() :
	string_key_(""),
	string_value_(""),
	real_value_(0.0),
	value_type_(string_value),
	ligand_chain_("")
{

}

AddJobPairData::AddJobPairData(AddJobPairData const & ) = default;

AddJobPairData::~AddJobPairData() = default;

void AddJobPairData::apply( Pose & pose)
{
	if ( ligand_chain_ == "" ) {
		if ( value_type_ == string_value ) {
			protocols::jd2::add_string_string_pair_to_current_job(string_key_,string_value_);
		} else if ( value_type_ == real_value ) {
			protocols::jd2::add_string_real_pair_to_current_job(string_key_,real_value_);
		} else {
			//This really shouldn't happen
			TR.Fatal << "No ligand found, and value_type is '" << value_type_ << "' for AddJobPairData" << std::endl;
			utility_exit_with_message("AddJobPairData needs either a ligand or a valid value type.");
		}
	} else {
		//Fundamental assumption being made that ligand is 1 residue.
		//This is true in the simple hts case but not universally true.
		core::Size chain_id = core::pose::get_chain_id_from_chain(ligand_chain_, pose);
		core::Size ligand_seqpos(pose.conformation().chain_begin(chain_id));
		core::chemical::ResidueType ligand_res_type(pose.conformation().residue_type(ligand_seqpos));
		if ( value_type_ == string_value ) {
			std::string value = ligand_res_type.get_string_property(string_key_);
			protocols::jd2::add_string_string_pair_to_current_job(string_key_, value);
		} else if ( value_type_ == real_value ) {
			core::Real value = ligand_res_type.get_numeric_property(string_key_);
			protocols::jd2::add_string_real_pair_to_current_job(string_key_,value);
		} else {
			//This really shouldn't happen
			TR.Fatal << "value_type is '" << value_type_ << "' for AddJobPairData" << std::endl;
			utility_exit_with_message("AddJobPairData needs a valid value type.");
		}
	}
}

// XRW TEMP std::string AddJobPairData::get_name() const
// XRW TEMP {
// XRW TEMP  return "AddJobPairData";
// XRW TEMP }

moves::MoverOP AddJobPairData::clone() const
{
	return moves::MoverOP( new AddJobPairData(*this) );
}

moves::MoverOP AddJobPairData::fresh_instance() const
{
	return moves::MoverOP( new AddJobPairData );
}

void AddJobPairData::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	if ( !tag->hasOption("value_type") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'value_type'");
	}
	if ( !tag->hasOption("key") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'key'");
	}
	if ( !tag->hasOption("value") && !tag->hasOption("value_from_ligand_chain") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'value' or 'value_from_ligand_chain'");
	}
	if ( tag->hasOption("value") && tag->hasOption("value_from_ligand_chain") ) {
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'value' or 'value_from_ligand_chain' but not both");
	}

	ligand_chain_ = tag->getOption<std::string>("value_from_ligand_chain","");

	std::string value_type_string(tag->getOption<std::string>("value_type"));

	if ( value_type_string == "string" ) {
		value_type_ = string_value;
	} else if ( value_type_string == "real" ) {
		value_type_ = real_value;
	} else {
		TR.Fatal << "Value type of '" << value_type_string << "' not supported in AddJobPairData." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' option 'value_type' can only be 'string' or 'real'");
	}

	string_key_ = tag->getOption<std::string>("key");
	if ( ligand_chain_ == "" ) {
		if ( value_type_ == string_value ) {
			string_value_ = tag->getOption<std::string>("value");
		} else if ( value_type_ == real_value ) {
			real_value_ = tag->getOption<core::Real>("value");
		} else {
			// If this happens all the error checking code above is wrong
			TR.Fatal << "Can't handle value type of '" << value_type_ << "'" << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' option 'value_type' can only be 'string' or 'real'");
		}
	}
}

std::string AddJobPairData::get_name() const {
	return mover_name();
}

std::string AddJobPairData::mover_name() {
	return "AddJobPairData";
}

void AddJobPairData::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//argument legality checker for "value_type":
	XMLSchemaRestriction value_type_type;
	value_type_type.name("value_type_type");
	value_type_type.base_type( xs_string );
	value_type_type.add_restriction( xsr_enumeration, "real");
	value_type_type.add_restriction( xsr_enumeration, "string");
	xsd.add_top_level_element( value_type_type );
	attlist + XMLSchemaAttribute::required_attribute( "value_type", "value_type_type", "type of value to add; must be 'string' or 'real'");

	attlist + XMLSchemaAttribute::required_attribute( "key", xs_string, "the string name for this item");

	std::string const value_warning("Of 'value' and 'value_from_ligand_chain', you must specify exactly one.");

	attlist + XMLSchemaAttribute( "value", xs_string, "Value to report; " + value_warning);
	attlist + XMLSchemaAttribute( "value_from_ligand_chain", xsct_char, "This should be a single letter describing the chain to get a value from.  It queries using the key." + value_warning);

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds either a String/Real or String/String pair to the JD2 Job for output", attlist );
}

std::string AddJobPairDataCreator::keyname() const {
	return AddJobPairData::mover_name();
}

protocols::moves::MoverOP
AddJobPairDataCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddJobPairData );
}

void AddJobPairDataCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddJobPairData::provide_xml_schema( xsd );
}


}
}
