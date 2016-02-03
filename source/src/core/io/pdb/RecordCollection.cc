// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/pdb/RecordCollection.cc
/// @brief   Method definitions for RecordCollection.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/pdb/RecordCollection.hh>
#include <core/io/pdb/record_def_io.hh>
#include <core/io/pdb/Field.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/string_util.hh>
#include <utility/exit.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>

#ifdef MULTI_THREADED
#ifdef CXX11

// C++11 headers
#include <atomic>
#include <mutex>

#endif
#endif


// Singleton set-up
namespace utility {

using core::io::pdb::RecordCollection;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< RecordCollection >::singleton_mutex_ {};
template <> std::atomic< RecordCollection * > utility::SingletonBase< RecordCollection >::instance_( 0 );
#else
template <> RecordCollection * utility::SingletonBase< RecordCollection >::instance_( 0 );
#endif

}  // namespace utility


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.io.pdb.RecordCollection" );


namespace core {
namespace io {
namespace pdb {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access
bool
RecordCollection::is_valid_record_type( std::string const & type )
{
	return get_instance()->get_record_definitions_map().count( utility::trim( type ) );
}

Record
RecordCollection::record_from_record_type( std::string const & type )
{
	std::string const & trimmed_type( utility::trim( type ) );
	if ( ! is_valid_record_type( trimmed_type ) ) {
		return get_instance()->get_record_definitions_map()[ "UNKNOW" ];
	}
	return get_instance()->get_record_definitions_map()[ trimmed_type ];
}

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RecordCollection::RecordCollection() :
	record_definitions_( read_record_definitions_from_file(
	basic::database::full_name( "input_output/pdb_record_defs" ) ) )
{}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
RecordCollection *
RecordCollection::create_singleton_instance()
{
	return new RecordCollection;
}

RecordRef
RecordCollection::get_record_definitions_map() const
{
	return record_definitions_;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
