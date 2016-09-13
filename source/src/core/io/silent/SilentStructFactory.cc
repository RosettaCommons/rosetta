// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/SilentStructFactory.cc
/// @brief Factory for creating various types of silent.
/// @author Greg Taylor <gktaylor@u.washington.edu>

// Unit headers
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/util.hh>

// Package headers
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructCreator.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace io {
namespace silent {

static THREAD_LOCAL basic::Tracer tr( "core.io.silent" );

/// @details Private constructor insures correctness of singleton.
SilentStructFactory::SilentStructFactory() {}

SilentStructFactory *
SilentStructFactory::create_singleton_instance()
{
	return new SilentStructFactory;
}

void
SilentStructFactory::factory_register( SilentStructCreatorCOP creator )
{
	ss_types_[ creator->keyname() ] = creator;
}

void
SilentStructFactory::show_available_silent_struct_types(
	std::ostream & out
) {
	out << "Available silent struct types:" << std::endl;
	for (
			std::map< std::string, io::silent::SilentStructCreatorCOP >::const_iterator
			it = ss_types_.begin(), ite = ss_types_.end(); it != ite; ++it
			) {
		out << "\t" << it->first << std::endl;
	}
}


bool
SilentStructFactory::has_silent_struct_type(
	std::string const & type_name
) {
	return ss_types_.find(type_name) != ss_types_.end();
}

SilentStructOP
SilentStructFactory::get_silent_struct( std::string const & type_name )
{
	tr.Trace << "generate silent struct of type " << type_name << std::endl;
	SilentStructCreatorMap::const_iterator iter = ss_types_.find( type_name );
	if ( iter != ss_types_.end() ) {
		return iter->second->create_silent_struct();
	} else {

		using std::string;
		using utility::vector1;
		string msg("SilentStructFactory::get_instance()->get_silent_struct:  ");
		msg += type_name + " does not name a known SilentStructType --> " +
			"check spelling or register new SilentStruct type in SilentStructFactory!";
		msg += "known types are:\n";
		vector1< string > types = get_ss_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}
		msg += "If you were using binary_rna, switch to binary\n";

		utility_exit_with_message( msg );
	}
	return 0;
}

utility::vector1< std::string >
SilentStructFactory::get_ss_names() const {
	using std::string;
	using utility::vector1;

	vector1< string > ss_names;
	for ( SilentStructCreatorMap::const_iterator
			it = ss_types_.begin(), end = ss_types_.end(); it != end; ++it ) {
		ss_names.push_back( it->first );
	}

	return ss_names;
}

SilentStructCreatorCOP
SilentStructFactory::get_creator( std::string const & type_name )
{
	SilentStructCreatorMap::const_iterator iter = ss_types_.find( type_name );
	if ( iter != ss_types_.end() ) {
		return iter->second;
	} else {

		using std::string;
		using utility::vector1;
		string msg("SilentStructFactory::get_creator:  ");
		msg += type_name + " does not name a known SilentStructType --> " +
			"check spelling or register new SilentStruct type in SilentStructFactory!";

		msg += "known types are:\n";
		vector1< string > types = get_ss_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}

		utility_exit_with_message( msg );
	}
	return 0;
}

SilentStructOP SilentStructFactory::get_silent_struct_in() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_silent_struct( option[ in::file::silent_struct_type ]() );
}

SilentStructOP SilentStructFactory::get_silent_struct_out() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	return get_silent_struct( option[ out::file::silent_struct_type ]() );
}

SilentStructOP SilentStructFactory::get_silent_struct_out( core::pose::Pose const& pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	bool already_binary( option[ out::file::silent_struct_type ]().find( "binary" ) != std::string::npos );
	bool const score_only (option[ out::file::silent_struct_type ]() == "score");
	bool const score_jump(option[ out::file::silent_struct_type ]() == "score_jump");

	// if the user has explicitly set the output type, honor the request
	if ( option[ out::file::silent_struct_type ].user() ) {
		return get_silent_struct( option[ out::file::silent_struct_type ]() ) ;
	}

	// if not, decide for the user or use the default
	if ( !score_only && !score_jump && !already_binary && !core::pose::is_ideal_pose(pose) ) {
		std::string binary_string( "binary" );
		if ( option[ out::file::silent_struct_type ]() == "rna" ) {
			binary_string = "binary_rna";
		};
		tr.Info << "detected attempt to write non-ideal pose to silent-file..."
			<< "Automatically switching to " << binary_string << " silent-struct type" << std::endl;
		return get_silent_struct( binary_string );
	}

	// use default
	return get_silent_struct( option[ out::file::silent_struct_type ]() ) ;
}


} // silent
} // io
} // core
