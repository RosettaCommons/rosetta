// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/splice/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarel@weizmann.ac.il)

// Unit headers
#include <devel/splice/SpliceSegment.hh>
#include <core/sequence/SequenceProfile.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
// Package headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <core/pose/util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <sstream>
#include <utility/string_util.hh>
#include <dirent.h>

namespace devel {
namespace splice {

static basic::Tracer TR( "devel.splice.SpliceSegment" );

using namespace core::sequence;
using namespace std;

SpliceSegment::SpliceSegment(){
	sequence_profile_.clear();
	pdb_to_profile_map_.clear();
}

SpliceSegment::~SpliceSegment() {}


void//rapper function for read_profile to accomedate the new way of reading profiles
SpliceSegment::read_many( string const Protein_family_path , string const segment){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string fname;
	for(size_t i = 1, i_end = option[ in::path::database ]().size(); i <= i_end; ++i) {
		fname = option[ in::path::database ](i).name();
		//TR<<"Trying to open folder:"<<fname<<std::endl;
	}


	//Get all PDB_profile_Match files from raltive path
	const std::string target_path = (fname+Protein_family_path+"pssm/"+segment+"/");


	DIR *dir;
	struct dirent *ent;
	const char * c =target_path.c_str();
	if ((dir = opendir (c))!= NULL) {
		while ((ent = readdir (dir)) != NULL) {
			//TR<<"Loading PSSM from file "<<target_path+ent->d_name<<std::endl;
			utility::io::izstream data( target_path+ent->d_name);
			if ( !data )
				continue;
			/*get segment name from file name*/
			string segment_name(target_path+ent->d_name);
			// Remove directory if present.
			// Do this before extension removal incase directory has a period character.
			const size_t last_slash_idx = segment_name.find_last_of("\\/");
			if (std::string::npos != last_slash_idx)
			{
				segment_name.erase(0, last_slash_idx + 1);
			}

			// Remove extension if present.
			const size_t period_idx = segment_name.rfind('.');
			if (std::string::npos != period_idx)
			{
				segment_name.erase(period_idx);
			}
		//	TR<<"segment name:"<<segment_name<<std::endl;
			if ((segment_name.compare("") == 0)||(segment_name.compare(".") == 0)){//hack to avoid trying to read . or.. when looking into directory
				continue;
			}
			read_profile(target_path+ent->d_name,segment_name);
		}
		closedir (dir);
	}

	else {
		/* could not open directory */
		utility_exit_with_message( "Directory path" + target_path+"was not found" );

	}
}

void //this version of read_profile is an older one and it allows the user to enter it's own PSSM files, instead of using the ones in Rosetta_DB
SpliceSegment::read_profile( string const file_name, string const segment_name ){
	SequenceProfileOP new_seqprof( new SequenceProfile );
	utility::file::FileName const fname( file_name );
	//TR<<"reading sequence profile from "<<file_name<<std::endl;
	new_seqprof->read_from_file( fname );
	//TR<< "SpliceSegmentsSeqProf:"<< new_seqprof->prof_row(1)<<std::endl;
	//TR<< "SpliceSegmentspprobabilty:"<< new_seqprof->probability_row(1)<<std::endl;
	sequence_profile_.insert( pair< string, SequenceProfileOP >( segment_name, new_seqprof ) );
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
	//TR<<"size of sequence_Profile_ is:"<<sequence_profile_.size()<<std::endl;
	//TR<<"pdb to profile map is:"<<pdb_name<<":"<<pdb_to_profile_map_[pdb_name]<<std::endl;
	return sequence_profile_[ pdb_to_profile_map_[ pdb_name ] ];
}

core::sequence::SequenceProfileOP
SpliceSegment::sequence_profile( string const profile_name ){
	return sequence_profile_[ profile_name ];
}

/// @brief generate one long sequence profile out of several sequence-profile segments.
SequenceProfileOP
concatenate_profiles( utility::vector1< SequenceProfileOP > const profiles, utility::vector1< std::string > segment_names_ordered){
	SequenceProfileOP concatenated_profile( new SequenceProfile );
	core::Size current_profile_size( 0 );
	core::Size current_segment_name( 1 );
	foreach( SequenceProfileOP const prof, profiles ){
		TR<<"now adding profile of segment "<< segment_names_ordered[ current_segment_name]<<std::endl;
		for( core::Size pos = 1; pos <= prof->size(); ++pos ){
			current_profile_size++;
			//TR<<"The sequence profile row for this residue is: "<<prof->prof_row(pos)<<std::endl;
		
			concatenated_profile->prof_row( prof->prof_row( pos ), current_profile_size );
			if (basic::options::option[ basic::options::OptionKeys::out::file::occurrence_data ].value()){
				//this option is by default false
				concatenated_profile->probabilty_row( prof->probability_row( pos ), current_profile_size );
			}

		}//for
		++current_segment_name;
	}//foreach
	TR<<"concatenated "<<profiles.size()<<" profiles with a total length of "<<current_profile_size<<std::endl;
	return concatenated_profile;
}

/// @brief reads the pdb-profile match file, matching each pdb file to one profile
void
SpliceSegment::read_pdb_profile( std::string const file_name ){
	utility::io::izstream data( file_name );
	if ( !data )
		utility_exit_with_message( "File not found " + file_name );
	//TR<<"Loading pdb profile pairs from file "<<file_name<<std::endl;
	string line;
	while (getline( data, line )){
		istringstream line_stream( line );
		while( !line_stream.eof() ){
			std::string pdb, profile;
			line_stream >> pdb >> profile;
			pdb_to_profile_map_.insert( pair< string, string >( pdb, profile ) );
			//TR<<"Loading pdb-profile pair: "<<pdb<<" "<<profile<<std::endl;
		}
	}

}

void
SpliceSegment::all_pdb_profile( string const Protein_family_path, string const segment){

	//This is to get the path to the database
	utility::io::izstream data;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string fname;
	for(size_t i = 1, i_end = option[ in::path::database ]().size(); i <= i_end; ++i) {
		fname = option[ in::path::database ](i).name();
		//TR<<"Trying to open folder:"<<fname<<std::endl;
	}


	//Get all PDB_profile_Match files from raltive path
	const std::string target_path = (fname+Protein_family_path+"pdb_profile_match/");


	DIR *dir;
	struct dirent *ent;
	const char * c =target_path.c_str();
	if ((dir = opendir (c))!= NULL) {
			string fileName(target_path+"pdb_profile_match."+segment);
			utility::io::izstream data( fileName );
				if ( !data )
					utility_exit_with_message( "File not found " + fileName );
			std::string line;
			TR<<"Loading pdb profile pairs from file "<<fileName<<std::endl;
			while (getline( data, line )){
				istringstream line_stream( line );
				while( !line_stream.eof() ){
					std::string pdb, profile;
					line_stream >> pdb >> profile;
					pdb_to_profile_map_.insert( pair< string, string >( pdb, profile ) );
					//TR<<"Loading pdb-profile pair: "<<pdb<<" "<<profile<<std::endl;
				}
			}

		closedir (dir);
	} else {
		/* could not open directory */
		utility_exit_with_message( "Directory path" + target_path+"was not found" );

	}
}


} //splice
} //devel
