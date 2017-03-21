// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/movers/AddorRemoveCsts.cc
/// @brief
/// @author Florian Richter (floric@u.washington.edu)

//unit headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/AddOrRemoveMatchCstsCreator.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>

//project headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/scoring/ScoreFunction.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// amw debug
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace enzdes {


//initializing static member variable
std::map< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP > AddOrRemoveMatchCsts::cstfile_map_;

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.AddorRemoveCsts" );

AddOrRemoveMatchCsts::AddOrRemoveMatchCsts()
: Mover("AddOrRemoveMatchCsts"),
	option_cstfile_(""), cstfile_(""),
	cst_action_(VOID), keep_covalent_(false),
	accept_blocks_missing_header_(false), fail_on_constraints_missing_(true)
{
	if ( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ) {
		option_cstfile_ = basic::options::option[basic::options::OptionKeys::enzdes::cstfile].value();
	}
	sfxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );

	//this might change in the future to a scorefunction that's fully definable
	//through a weights file, but atm it's probably safe to simply set the constraint
	//weights to non zero values
	sfxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfxn_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfxn_->set_weight( core::scoring::angle_constraint, 1.0 );
	sfxn_->set_weight( core::scoring::dihedral_constraint, 1.0 );
}

AddOrRemoveMatchCsts::AddOrRemoveMatchCsts( AddOrRemoveMatchCsts const & ) = default;

AddOrRemoveMatchCsts::~AddOrRemoveMatchCsts()= default;

protocols::moves::MoverOP
AddOrRemoveMatchCsts::clone() const
{
	return protocols::moves::MoverOP( new AddOrRemoveMatchCsts( *this ) );
}

protocols::moves::MoverOP
AddOrRemoveMatchCsts::fresh_instance() const
{
	return protocols::moves::MoverOP( new AddOrRemoveMatchCsts() );
}

void
AddOrRemoveMatchCsts::apply( core::pose::Pose & pose )
{
	std::string cstfile( cstfile_ != "" ? cstfile_ : option_cstfile_ );
	//safety check for this function being called without a cstfile having been specified
	if ( (cstfile == "") && (cstfile_map_.size() == 0) ) {
		tr.Warning << "apply function of enzdes constraints mover called even though no cstfile has been specified on the commandline, in the tag, or programmatically. This function will have no effect." << std::endl;
		return;
	}

	toolbox::match_enzdes_util::EnzConstraintIOOP cst_io( this->get_EnzConstraintIO_for_cstfile( cstfile ) );

	switch (cst_action_)  {
	case ADD_NEW :
		cst_io->add_constraints_to_pose( pose, sfxn_, accept_blocks_missing_header_ );
		break;
	case ADD_PREGENERATED :
		cst_io->add_pregenerated_constraints_to_pose( pose, sfxn_ );
		break;
	case REMOVE :
		cst_io->remove_constraints_from_pose( pose, keep_covalent_, fail_on_constraints_missing_ );
		break;
	case VOID :
		utility_exit_with_message("Illegal action for AddOrRemoveMatchCsts mover specified.");
		break;
	}
}

void
AddOrRemoveMatchCsts::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	cstfile_ = tag->getOption<std::string>( "cstfile", "" );
	if ( (cstfile_ == "") && ( option_cstfile_ == "" ) ) {
		tr.Warning << "No name for the enzdes .cst file was specified in either the options, the xml tag, or programatically. AddOrRemoveMatchCsts will turn into a null operation." << std::endl;
	}

	std::string cst_instruction = tag->getOption<std::string>( "cst_instruction", "void" );
	if ( cst_instruction == "add_new" ) cst_action_ = ADD_NEW;
	else if ( cst_instruction == "add_pregenerated" ) cst_action_ = ADD_PREGENERATED;
	else if ( cst_instruction == "remove" ) cst_action_ = REMOVE;
	else {
		throw utility::excn::EXCN_RosettaScriptsOption("Illegal or no value for cst_instruction in xml tag given. Has to be either 'add_new', 'add_pregenerated', or 'remove'.");
	}

	keep_covalent_ = tag->getOption<bool>( "keep_covalent", 0 );
	accept_blocks_missing_header_ = tag->getOption<bool>( "accept_blocks_missing_header", 0 );
	fail_on_constraints_missing_ = tag->getOption<bool>( "fail_on_constraints_missing", 1 );
}

void
AddOrRemoveMatchCsts::cstfile( std::string const & setting )
{
	cstfile_ = setting;
}


toolbox::match_enzdes_util::EnzConstraintIOCOP
AddOrRemoveMatchCsts::get_const_EnzConstraintIO_for_cstfile( std::string cstfile )
{

	if ( cstfile == "" ) {
		if ( cstfile_map_.size() == 1 ) return cstfile_map_.begin()->second;
	}
	std::map< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP>::const_iterator enzio_it = cstfile_map_.find( cstfile );
	if ( enzio_it == cstfile_map_.end() ) {
		tr.Warning << "trying to get an EnzConstraintIOOP object for cstfile " << cstfile << " that hasn't been instantiated yet. Returning NULL pointer." << std::endl;
		return nullptr;
	}
	return enzio_it->second;

}

toolbox::match_enzdes_util::EnzConstraintIOOP
AddOrRemoveMatchCsts::get_EnzConstraintIO_for_cstfile(
	std::string const & cstfile )
{
	auto cstio_it = cstfile_map_.find( cstfile );

	if ( cstio_it == cstfile_map_.end() ) {
		toolbox::match_enzdes_util::EnzConstraintIOOP new_cst_io( new toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
		new_cst_io->read_enzyme_cstfile( cstfile );
		cstfile_map_.insert( std::pair< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP >( cstfile, new_cst_io ) );
		return new_cst_io;
	}
	return cstio_it->second;
}

std::string AddOrRemoveMatchCsts::get_name() const {
	return mover_name();
}

std::string AddOrRemoveMatchCsts::mover_name() {
	return "AddOrRemoveMatchCsts";
}

void AddOrRemoveMatchCsts::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"cstfile", xs_string,
		"name of file to get csts from (can be specified here if one wants to change the constraints, "
		"e.g. tighten or relax them, as the pose progresses down a protocol.)");

	XMLSchemaRestriction cst_instr_enum;
	cst_instr_enum.name("cst_instruction_types");
	cst_instr_enum.base_type(xs_string);
	cst_instr_enum.add_restriction(xsr_enumeration, "add_new");
	cst_instr_enum.add_restriction(xsr_enumeration, "add_pregenerated");
	cst_instr_enum.add_restriction(xsr_enumeration, "remove");
	xsd.add_top_level_element( cst_instr_enum );

	attlist + XMLSchemaAttribute::required_attribute(
		"cst_instruction", "cst_instruction_types",
		"1 of 3 choices - \"add_new\" (read from file), \"remove\", or \"add_pregenerated\" "
		"(i.e. if enz csts existed at any point previosuly in the protocol add them back)");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"keep_covalent", xsct_rosetta_bool,
		"during removal, keep constraints corresponding to covalent bonds between protein and ligand intact (default=0).",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"accept_blocks_missing_header", xsct_rosetta_bool,
		"allow more blocks in the cstfile than specified in header REMARKs (see enzdes documentation for details, default=0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"fail_on_constraints_missing", xsct_rosetta_bool,
		"When removing constraints, raise an error if the constraint blocks do not exist in the pose (default=1).",
		"true");


	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"XRW XSD: TO DO",
		attlist );
}

std::string AddOrRemoveMatchCstsCreator::keyname() const {
	return AddOrRemoveMatchCsts::mover_name();
}

protocols::moves::MoverOP
AddOrRemoveMatchCstsCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddOrRemoveMatchCsts );
}

void AddOrRemoveMatchCstsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddOrRemoveMatchCsts::provide_xml_schema( xsd );
}


} //enzdes
} //protocols
