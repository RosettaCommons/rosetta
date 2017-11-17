// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/LoopsFileFallbackConfiguration.cc
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit Headers
#include <protocols/loops/LoopsFileFallbackConfiguration.hh>
#include <protocols/loops/LoopsFileFallbackConfigurationCreator.hh>


// Platform Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// basic headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// numeric headers
#include <numeric/random/random.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>

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


LoopsFileFallbackConfiguration::LoopsFileFallbackConfiguration()
{}

/// @details Return true if the user has set the "-loops:loop_file" flag and false otherwise.
bool
LoopsFileFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user();
}

/// @details The return value is "loops_file" to indicate that a LoopsFileLoader is required.
basic::resource_manager::LoaderType
LoopsFileFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "LoopsFile";
}

/// @details The %locator_id for the fallback configuration is set by the options system.  The
/// get_loops_filename_from_options() helper method is used to handle complex cases.
basic::resource_manager::LocatorID
LoopsFileFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return get_loops_filename_from_options();
}

/// @details Return a NULL pointer to trigger the creation of a default LoopsFileOptions later in the %resource creation
/// process.
basic::resource_manager::ResourceOptionsOP
LoopsFileFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// use the default loops file options.
	return nullptr;
}

/// @details Return a string that provides a helpful message to the user so s/he can determine how to correctly use the
/// loops_file resource.
std::string
LoopsFileFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The LoopsFileFallbackConfiguration requires that the flag '-loops:loop_file' be set on the command line.";
}

/// @details There are three options system scenarios this method can handle.  Each scenario and the corresponding
///behavior is outlined below:
/// @li A single filename is specified - return the filename as a string
/// @li Several filenames are specified - return one of the filenames, chosen at random, as a string
/// @li No filename is specified - throw an exception requesting that the user double check his/her flags.
/// @throws EXCN_Msg_Exception
basic::resource_manager::LocatorID
LoopsFileFallbackConfiguration::get_loops_filename_from_options() const
{
	// the next line uses value_or to avoid a call to std::exit in the options class.  I can test things that throw exceptions.  Just sayin'.
	utility::vector1< std::string > loops_files = basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value_or( utility::vector1< std::string >() );
	if ( ! loops_files.size() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "The fallback LoopsFile resource option has no loops files associated with it! Was the option omitted from the command line?");
	}
	core::Size const which_loops_file( loops_files.size() == 1 ? 1 : core::Size( numeric::random::rg().random_range(1,( loops_files.size() ))));
	return loops_files[ which_loops_file ];
}

/// @details Return an owning pointer to a newly constructed default instance of FallbackConfiguration.
basic::resource_manager::FallbackConfigurationOP
LoopsFileFallbackConfigurationCreator::create_fallback_configuration() const
{
	return basic::resource_manager::FallbackConfigurationOP( new LoopsFileFallbackConfiguration );
}

/// @details Return a string specifying the type of %FallbackConfiguration to create (loops_file).
std::string
LoopsFileFallbackConfigurationCreator::resource_description() const
{
	return "LoopsFile";
}

} // namespace loops
} // namespace protocols
