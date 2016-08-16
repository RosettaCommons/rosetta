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

#include <protocols/simple_moves/AddJobPairData.hh>
#include <protocols/simple_moves/AddJobPairDataCreator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace simple_moves {

std::string AddJobPairDataCreator::keyname() const
{
	return AddJobPairDataCreator::mover_name();
}

moves::MoverOP AddJobPairDataCreator::create_mover() const
{
	return moves::MoverOP( new AddJobPairData );
}

std::string AddJobPairDataCreator::mover_name()
{
	return "AddJobPairData";
}

AddJobPairData::AddJobPairData() :
	string_key_(""),
	string_value_(""),
	real_value_(0.0),
	value_type_(string_value),
	ligand_chain_("")
{

}

AddJobPairData::AddJobPairData(AddJobPairData const & src) : Mover(src),
	string_key_(src.string_key_),
	string_value_(src.string_value_),
	real_value_(src.real_value_),
	value_type_(src.value_type_),
	ligand_chain_(src.ligand_chain_)
{

}

AddJobPairData::~AddJobPairData()
{

}

void AddJobPairData::apply( Pose & pose)
{
	jd2::JobOP current_job = jd2::JobDistributor::get_instance()->current_job();

	if ( ligand_chain_ == "" ) {
		if ( value_type_ == string_value ) {
			current_job->add_string_string_pair(string_key_,string_value_);
		} else if ( value_type_ == real_value ) {
			current_job->add_string_real_pair(string_key_,real_value_);
		} else {
			//This really shouldn't happen
			assert(false);
		}
	} else {
		//Fundamental assumption being made that ligand is 1 residue.
		//This is true in the simple hts case but not universally true.
		core::Size chain_id = core::pose::get_chain_id_from_chain(ligand_chain_, pose);
		core::Size ligand_seqpos(pose.conformation().chain_begin(chain_id));
		core::chemical::ResidueType ligand_res_type(pose.conformation().residue_type(ligand_seqpos));
		if ( value_type_ == string_value ) {
			std::string value = ligand_res_type.get_string_property(string_key_);
			current_job->add_string_string_pair(string_key_, value);
		} else if ( value_type_ == real_value ) {
			core::Real value = ligand_res_type.get_numeric_property(string_key_);
			current_job->add_string_real_pair(string_key_,value);
		} else {
			//This really shouldn't happen
			assert(false);
		}
	}
}

std::string AddJobPairData::get_name() const
{
	return "AddJobPairData";
}

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
			assert(false);
		}
	}
}

}
}
