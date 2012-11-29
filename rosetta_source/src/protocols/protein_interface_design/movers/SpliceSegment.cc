// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarel@weizmann.ac.il)

// Unit headers
#include <protocols/protein_interface_design/movers/SpliceSegment.hh>
#include <core/sequence/SequenceProfile.hh>
#include <utility/exit.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
// Package headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <sstream>

namespace protocols {
namespace protein_interface_design {
namespace movers {

static basic::Tracer TR( "protocols.protein_interface_design.movers.SpliceSegment" );

using namespace core::sequence;
using namespace std;

SpliceSegment::SpliceSegment(){
	sequence_profile_.clear();
	pdb_to_profile_map_.clear();
}

SpliceSegment::~SpliceSegment() {}

void
SpliceSegment::read_profile( string const file_name, string const segment_name ){
	SequenceProfileOP new_seqprof( new SequenceProfile );
	utility::file::FileName const fname( file_name );
	TR<<"reading sequence profile from "<<file_name<<std::endl;
	new_seqprof->read_from_file( fname );
	sequence_profile_.insert( pair< string, SequenceProfileOP >( file_name, new_seqprof ) );
}

SequenceProfileOP
SpliceSegment::get_profile( std::string const segment_name ){
	return sequence_profile_[ segment_name ];
}

void
SpliceSegment::add_pdb_profile_pair( std::string const pdb, string const profile_name ){
	pdb_to_profile_map_.insert( pair< string, string >( pdb, profile_name ) );
}

/// @brief this is the most useful method, returning the sequence profile according to a pdb file name
SequenceProfileOP
SpliceSegment::pdb_profile( std::string const pdb_name ){
	return sequence_profile_[ pdb_to_profile_map_[ pdb_name ] ];
}

core::sequence::SequenceProfileOP
SpliceSegment::sequence_profile( string const profile_name ){
	return sequence_profile_[ profile_name ];
}

/// @brief generate one long sequence profile out of several sequence-profile segments.
SequenceProfileOP
concatenate_profiles( utility::vector1< SequenceProfileOP > const profiles ){
	SequenceProfileOP concatenated_profile( new SequenceProfile );
	core::Size current_profile_size( 0 );
	foreach( SequenceProfileOP const prof, profiles ){
		for( core::Size pos = 1; pos <= prof->size(); ++pos ){
			current_profile_size++;
			concatenated_profile->prof_row( prof->prof_row( pos ), current_profile_size );
		}
	}
	TR<<"concatenated "<<profiles.size()<<" profiles with a total length of "<<concatenated_profile->size()<<std::endl;
}

/// @brief reads the pdb-profile match file, matching each pdb file to one profile
void
SpliceSegment::read_pdb_profile( std::string const file_name ){
  utility::io::izstream data( file_name );
  if ( !data )
		utility_exit_with_message( "File not found " + file_name );
	TR<<"Loading pdb profile pairs from file "<<file_name<<std::endl;
	string line;
	getline( data, line );
	if( line.length() > 2 )
		utility_exit_with_message( "pdb/profile file empty or corrupted "+file_name );
	istringstream line_stream( line );
	while( !line_stream.eof() ){
		std::string pdb, profile;
		line_stream >> pdb >> profile;
		pdb_to_profile_map_.insert( pair< string, string >( pdb, profile ) );
	}
}

} //movers
} //protein_interface_design
} //protocols
