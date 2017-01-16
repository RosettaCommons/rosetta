// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileRMOptions.hh
/// @brief  Options for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/io/silent/SilentFileRMOptionsCreator.hh>
#include <core/io/silent/SilentFileRMOptions.hh>
#include <utility/tag/Tag.hh>

// Project Headers
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileOptions.hh>

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

///// SilentFileRMOptionsCreator /////
SilentFileRMOptionsCreator::SilentFileRMOptionsCreator() {}

SilentFileRMOptionsCreator::~SilentFileRMOptionsCreator() {}

ResourceOptionsOP
SilentFileRMOptionsCreator::create_options() const {
	return ResourceOptionsOP( new SilentFileRMOptions );
}

string
SilentFileRMOptionsCreator::options_type() const {
	return SilentFileRMOptions::class_name();
}

SilentFileRMOptions::SilentFileRMOptions() :
	ResourceOptions(),
	options_( new SilentFileOptions )
{}

SilentFileRMOptions::SilentFileRMOptions(
	string const & name
) :
	ResourceOptions(name),
	options_( new SilentFileOptions )
{}

SilentFileRMOptions::~SilentFileRMOptions() {}

SilentFileRMOptions::SilentFileRMOptions(
	SilentFileRMOptions const & src
) :
	ResourceOptions(src),
	options_( new SilentFileOptions( *src.options_) )
{}

SilentFileOptions const &
SilentFileRMOptions:: opts() const
{
	return *options_;
}

void
SilentFileRMOptions::parse_my_tag(
	TagCOP tag
)
{
	options_->read_from_tag( tag );
}

std::string
SilentFileRMOptions::type() const { return class_name(); }

std::string
SilentFileRMOptions::class_name()  { return "SilentFileOptions"; }

} // namespace
} // namespace
} // namespace
