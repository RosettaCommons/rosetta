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


// Unit headers
#include <core/io/pdb/RecordCollection.hh>
#include <core/io/pdb/record_def_io.hh>
#include <core/io/pdb/Field.hh>

// Project header
#include <core/types.hh>

// Utility headers
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
	return get_instance()->string_to_record_type_map_.count( utility::trim( type ) );
}

Record
RecordCollection::record_from_record_type( RecordType const & type )
{
	return get_instance()->record_definitions_[ type ];
}

Record
RecordCollection::record_from_record_type( std::string const & type )
{
	RecordCollection * instance( get_instance() );
	if ( ! is_valid_record_type( type ) ) {
		return instance->record_definitions_[ UNKNOW ];
	}
	return instance->record_definitions_[ instance->string_to_record_type_map_[ utility::trim( type ) ] ];
}

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RecordCollection::RecordCollection()
{
	string_to_record_type_map_[ "ANISOU" ] = ANISOU;
	string_to_record_type_map_[ "ATOM" ] = ATOM;
	string_to_record_type_map_[ "AUTHOR" ] = AUTHOR;
	string_to_record_type_map_[ "CAVEAT" ] = CAVEAT;
	string_to_record_type_map_[ "CISPEP" ] = CISPEP;
	string_to_record_type_map_[ "COMPND" ] = COMPND;
	string_to_record_type_map_[ "CONECT" ] = CONECT;
	string_to_record_type_map_[ "CRYST1" ] = CRYST1;
	string_to_record_type_map_[ "DBREF" ] = DBREF;
	string_to_record_type_map_[ "DBREF1" ] = DBREF1;
	string_to_record_type_map_[ "DBREF2" ] = DBREF2;
	string_to_record_type_map_[ "END" ] = END;
	string_to_record_type_map_[ "ENDMDL" ] = ENDMDL;
	string_to_record_type_map_[ "EXPDTA" ] = EXPDTA;
	string_to_record_type_map_[ "FORMUL" ] = FORMUL;
	string_to_record_type_map_[ "HEADER" ] = HEADER;
	string_to_record_type_map_[ "HELIX" ] = HELIX;
	string_to_record_type_map_[ "HET" ] = HET;
	string_to_record_type_map_[ "HETATM" ] = HETATM;
	string_to_record_type_map_[ "HETNAM" ] = HETNAM;
	string_to_record_type_map_[ "HETSYN" ] = HETSYN;
	string_to_record_type_map_[ "JRNL" ] = JRNL;
	string_to_record_type_map_[ "KEYWDS" ] = KEYWDS;
	string_to_record_type_map_[ "LINK" ] = LINK;
	string_to_record_type_map_[ "MASTER" ] = MASTER;
	string_to_record_type_map_[ "MDLTYP" ] = MDLTYP;
	string_to_record_type_map_[ "MODEL" ] = MODEL;
	string_to_record_type_map_[ "MODRES" ] = MODRES;
	string_to_record_type_map_[ "MTRIX1" ] = MTRIX1;
	string_to_record_type_map_[ "MTRIX2" ] = MTRIX2;
	string_to_record_type_map_[ "MTRIX3" ] = MTRIX3;
	string_to_record_type_map_[ "NUMMD" ] = NUMMD;
	string_to_record_type_map_[ "OBSLTE" ] = OBSLTE;
	string_to_record_type_map_[ "ORIGX1" ] = ORIGX1;
	string_to_record_type_map_[ "ORIGX2" ] = ORIGX2;
	string_to_record_type_map_[ "ORIGX3" ] = ORIGX3;
	string_to_record_type_map_[ "REMARK" ] = REMARK;
	string_to_record_type_map_[ "REVDAT" ] = REVDAT;
	string_to_record_type_map_[ "SCALE1" ] = SCALE1;
	string_to_record_type_map_[ "SCALE2" ] = SCALE2;
	string_to_record_type_map_[ "SCALE3" ] = SCALE3;
	string_to_record_type_map_[ "SEQADV" ] = SEQADV;
	string_to_record_type_map_[ "SEQRES" ] = SEQRES;
	string_to_record_type_map_[ "SHEET" ] = SHEET;
	string_to_record_type_map_[ "SITE" ] = SITE;
	string_to_record_type_map_[ "SOURCE" ] = SOURCE;
	string_to_record_type_map_[ "SPLIT" ] = SPLIT;
	string_to_record_type_map_[ "SPRSDE" ] = SPRSDE;
	string_to_record_type_map_[ "SSBOND" ] = SSBOND;
	string_to_record_type_map_[ "TER" ] = TER;
	string_to_record_type_map_[ "TITLE" ] = TITLE;
	string_to_record_type_map_[ "UNKNOW" ] = UNKNOW;

	record_definitions_ = read_record_definitions_from_file(
		basic::database::full_name( "input_output/pdb_record_defs" ), string_to_record_type_map_ );
}

// Singleton-creation function for use with utility::thread::threadsafe_singleton
RecordCollection *
RecordCollection::create_singleton_instance()
{
	return new RecordCollection;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
