// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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


namespace protocols {
namespace enzdes {


//initializing static member variable
std::map< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP > AddOrRemoveMatchCsts::cstfile_map_;

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.AddorRemoveCsts" );

std::string
AddOrRemoveMatchCstsCreator::keyname() const
{
	return AddOrRemoveMatchCstsCreator::mover_name();
}

protocols::moves::MoverOP
AddOrRemoveMatchCstsCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddOrRemoveMatchCsts );
}

std::string
AddOrRemoveMatchCstsCreator::mover_name()
{
	return "AddOrRemoveMatchCsts";
}

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

AddOrRemoveMatchCsts::AddOrRemoveMatchCsts( AddOrRemoveMatchCsts const & other ) :
	Mover( other ),
	option_cstfile_(other.option_cstfile_), cstfile_(other.cstfile_),
	cst_action_(other.cst_action_), keep_covalent_(other.keep_covalent_),
	accept_blocks_missing_header_(other.accept_blocks_missing_header_),
	fail_on_constraints_missing_(other.fail_on_constraints_missing_),
	sfxn_(other.sfxn_)
{}

AddOrRemoveMatchCsts::~AddOrRemoveMatchCsts(){}

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
		tr.Warning << "Warning: apply function of enzdes constraints mover called even though no cstfile has been specified on the commandline, in the tag, or programmatically. This function will have no effect." << std::endl;
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

std::string
AddOrRemoveMatchCsts::get_name() const {
	return AddOrRemoveMatchCstsCreator::mover_name();
}

void
AddOrRemoveMatchCsts::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	cstfile_ = tag->getOption<std::string>( "cstfile", "" );
	if ( (cstfile_ == "") && ( option_cstfile_ == "" ) ) {
		tr.Warning << "WARNING: No name for the enzdes .cst file was specified in either the options, the xml tag, or programatically. AddOrRemoveMatchCsts will turn into a null operation." << std::endl;
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
		tr.Warning << "WARNING: trying to get an EnzConstraintIOOP object for cstfile " << cstfile << " that hasn't been instantiated yet. Returning NULL pointer." << std::endl;
		return NULL;
	}
	return enzio_it->second;

}

toolbox::match_enzdes_util::EnzConstraintIOOP
AddOrRemoveMatchCsts::get_EnzConstraintIO_for_cstfile(
	std::string const cstfile )
{
	std::map< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP >::iterator cstio_it = cstfile_map_.find( cstfile );

	if ( cstio_it == cstfile_map_.end() ) {
		toolbox::match_enzdes_util::EnzConstraintIOOP new_cst_io( new toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
		new_cst_io->read_enzyme_cstfile( cstfile );
		cstfile_map_.insert( std::pair< std::string, toolbox::match_enzdes_util::EnzConstraintIOOP >( cstfile, new_cst_io ) );
		return new_cst_io;
	}
	return cstio_it->second;
}

} //enzdes
} //protocols
