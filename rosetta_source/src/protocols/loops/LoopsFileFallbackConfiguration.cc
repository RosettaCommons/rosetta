// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/LoopsFileFallbackConfiguration.cc
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit Headers
#include <protocols/loops/LoopsFileFallbackConfiguration.hh>
#include <protocols/loops/LoopsFileFallbackConfigurationCreator.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// Platform Headers
#include <core/types.hh>
#include <utility/vector1.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>

// numeric headers
#include <numeric/random/random.hh>


//C++ Headers
#include <string>
#include <map>


namespace protocols {
namespace loops {

using basic::resource_manager::LoaderType;
using basic::resource_manager::LocatorID;
using basic::resource_manager::LocatorTag;
using basic::resource_manager::ResourceDescription;
using basic::resource_manager::ResourceTag;
using basic::resource_manager::ResourceOptionsTag;

static numeric::random::RandomGenerator RG(1337);

LoaderType LoopsFileFallbackConfiguration::loader_type_;
LocatorID LoopsFileFallbackConfiguration::locator_id_;
LocatorTag LoopsFileFallbackConfiguration::locator_tag_;
ResourceTag LoopsFileFallbackConfiguration::resource_tag_;
ResourceOptionsTag LoopsFileFallbackConfiguration::resource_options_tag_;


LoopsFileFallbackConfiguration::LoopsFileFallbackConfiguration()
{
	loader_type_ = "LoopsFile";
	locator_id_ = get_loops_filename_from_options();
	locator_tag_ = "";
	resource_tag_ = "LoopsFile_FallbackConfiguration";
	resource_options_tag_ = "";

}

ResourceTag const &
LoopsFileFallbackConfiguration::get_resource_tag_from_description( ResourceDescription const & ) const
{
	return resource_tag_;
}

LocatorTag const &
LoopsFileFallbackConfiguration::get_locator_tag_from_description( ResourceDescription const & ) const
{
	return locator_tag_;
}

LocatorID const &
LoopsFileFallbackConfiguration::get_locator_id_from_description( ResourceDescription const &) const
{
	return locator_id_;
}

LoaderType const &
LoopsFileFallbackConfiguration::get_loader_type_from_description( ResourceDescription const & ) const
{
	return loader_type_;
}

ResourceOptionsTag const &
LoopsFileFallbackConfiguration::get_resource_options_tag_from_description( ResourceDescription const & ) const
{
	return resource_options_tag_;
}
	

basic::resource_manager::LocatorID
LoopsFileFallbackConfiguration::get_loops_filename_from_options() const
{
	// the next line uses value_or to avoid a call to std::exit in the options class.  I can test things that throw exceptions.  Just sayin'.
	utility::vector1< std::string > loops_files = basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value_or( utility::vector1< std::string >() );
	if ( ! loops_files.size() )
	{
		throw utility::excn::EXCN_Msg_Exception("The fallback LoopsFile resource option has no loops files associated with it! Was the option omitted from the command line?");
	}
	core::Size const which_loops_file( loops_files.size() == 1 ? 1 : core::Size( RG.random_range(1,( loops_files.size() ))));
	return loops_files[ which_loops_file ];
}

basic::resource_manager::FallbackConfigurationOP
LoopsFileFallbackConfigurationCreator::create_fallback_configuration() const
{
	return new LoopsFileFallbackConfiguration;
}

std::string
LoopsFileFallbackConfigurationCreator::resource_description() const
{
	return "LoopsFile";
}

} // namespace loops
} // namespace protocols
