// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/sequence/SequenceFactory.cc
/// @brief Factory for creating various types of sequences.
/// @author James Thompson

#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceCreator.hh>
#include <core/sequence/SequenceFactory.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace sequence {

static THREAD_LOCAL basic::Tracer tr( "core.sequence.SequenceFactory" );

SequenceFactory *
SequenceFactory::create_singleton_instance()
{
	return new SequenceFactory;
}


/// @details Private constructor insures correctness of singleton.
SequenceFactory::SequenceFactory() {}

void
SequenceFactory::factory_register( SequenceCreatorCOP creator )
{
	seq_types_[ creator->keyname() ] = creator;
}

SequenceOP SequenceFactory::seq_from_file(
	std::string const & fn,
	std::string const & type_name
) {
	SequenceOP seq = get_sequence(type_name);
	tr.Debug << "reading from filename " << fn << std::endl;
	if ( utility::file::file_exists(fn) ) {
		seq->read_from_file(fn);
	} else {
		tr.Error << "Error: file " << fn << " doesn't exist!" << std::endl;
	}
	return seq;
}

SequenceOP
SequenceFactory::get_sequence( std::string const & type_name )
{
	tr.Trace << "generate sequence of type " << type_name << std::endl;
	SequenceCreatorMap::const_iterator iter = seq_types_.find( type_name );
	if ( iter != seq_types_.end() ) {
		return iter->second->create_sequence();
	} else {

		using std::string;
		using utility::vector1;
		string msg("SequenceFactory::get_instance()->get_sequence:  ");
		msg += type_name + " does not name a known SequenceType --> " +
			"check spelling or register new Sequence type in SequenceFactory!";

		msg += "known types are:\n";
		vector1< string > types = get_seq_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}

		utility_exit_with_message( msg );
	}
	return nullptr;
}

/// @details DANGER DANGER DANGER: NOT THREADSAFE!
/// This could be make threadsafe so long as seq_types_
/// is properly locked.
utility::vector1< std::string >
SequenceFactory::get_seq_names() const {
	using std::string;
	using utility::vector1;

	vector1< string > seq_names;
	for ( auto const & seq_type : seq_types_ ) {
		seq_names.push_back( seq_type.first );
	}

	return seq_names;
}

SequenceCreatorCOP
SequenceFactory::get_creator( std::string const & type_name )
{
	SequenceCreatorMap::const_iterator iter = seq_types_.find( type_name );
	if ( iter != seq_types_.end() ) {
		return iter->second;
	} else {

		using std::string;
		using utility::vector1;
		string msg("SequenceFactory::get_creator:  ");
		msg += type_name + " does not name a known SequenceType --> " +
			"check spelling or register new Sequence type in SequenceFactory!";

		msg += "known types are:\n";
		vector1< string > types = get_seq_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}

		utility_exit_with_message( msg );
	}
	return nullptr;
}

} // sequence
} // core
