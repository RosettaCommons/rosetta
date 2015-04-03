// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/silent/SilentFileOptions.hh
/// @brief  Options for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/io/silent/SilentFileOptionsCreator.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/io/silent/SilentStructFactory.hh>

// Plaform Headers
#include <core/types.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <string>
#include <sstream>

namespace core {
namespace io {
namespace silent {

using basic::resource_manager::ResourceOptionsOP;
using basic::resource_manager::ResourceOptions;
using std::endl;
using std::string;
using std::stringstream;
using utility::tag::TagCOP;

///// SilentFileOptionsCreator /////
SilentFileOptionsCreator::SilentFileOptionsCreator() {}

SilentFileOptionsCreator::~SilentFileOptionsCreator() {}

ResourceOptionsOP
SilentFileOptionsCreator::create_options() const {
	return ResourceOptionsOP( new SilentFileOptions );
}

string
SilentFileOptionsCreator::options_type() const {
	return "SilentFileOptions";
}

SilentFileOptions::SilentFileOptions() :
	ResourceOptions(),
	silent_struct_type_()
{}

SilentFileOptions::SilentFileOptions(
	string const & name
) :
	ResourceOptions(name),
	silent_struct_type_()
{}

SilentFileOptions::~SilentFileOptions() {}

SilentFileOptions::SilentFileOptions(
	SilentFileOptions const & src
) :
	ResourceOptions(src),
	silent_struct_type_(src.silent_struct_type_)
{}

std::string
SilentFileOptions::get_silent_struct_type() const {
	return silent_struct_type_;
}

void
SilentFileOptions::set_silent_struct_type(
	string const & silent_struct_type
) {
	if(!SilentStructFactory::get_instance()->has_silent_struct_type(
			silent_struct_type)){
		stringstream err_msg;
		err_msg
			<< "In SilentFileOptions,"
			<< " the silent_struct_type '" << silent_struct_type << "'"
			<< " was not recognized." << endl;
		SilentStructFactory::get_instance()->show_available_silent_struct_types(
			err_msg);
		throw utility::excn::EXCN_BadInput(err_msg.str());
	}

	silent_struct_type_ = silent_struct_type;
}

void
SilentFileOptions::parse_my_tag(
	TagCOP tag
) {
	if(!tag->hasOption("silent_struct_type")){
		throw utility::excn::EXCN_BadInput(
			"The SilentFileOptions requires the tag 'silent_struct_type',"
			" but it was not provided.");
	}

	silent_struct_type_ = tag->getOption<std::string>("silent_struct_type");

	if(!SilentStructFactory::get_instance()->has_silent_struct_type(
			silent_struct_type_)){
		stringstream err_msg;
		err_msg
			<< "In the SilentFileOptions with tag name '" << tag->getName() << "',"
			<< " the silent_struct_type '" << silent_struct_type_ << "'"
			<< " was not recognized." << endl;
		SilentStructFactory::get_instance()->show_available_silent_struct_types(
			err_msg);
		throw utility::excn::EXCN_BadInput(err_msg.str());
	}
}


} // namespace
} // namespace
} // namespace
