// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/splice/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarel@weizmann.ac.il)

// Unit headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <core/sequence/SequenceProfile.hh>
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
#include <basic/database/open.hh>

#ifdef WIN32
#include <dirent_windows.h>
#else
#include <dirent.h>
#endif

namespace protocols {
namespace splice {

static basic::Tracer TR( "protocols.splice.SpliceSegment" );

using namespace core::sequence;
using namespace std;


SpliceSegment::SpliceSegment()
: utility::pointer::ReferenceCount()
{
	sequence_profile_.clear();
	pdb_to_profile_map_.clear();
}

SpliceSegment::~SpliceSegment() = default;


void//Read in multiple pssm files
SpliceSegment::read_many( string const & Protein_family_path , string const & segment){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Get all PDB_profile_Match files from raltive path
	const std::string target_path = basic::database::full_name( Protein_family_path+"pssm/"+segment+"/" );


	DIR *dir;
	const char * c =target_path.c_str();
	if ( (dir = opendir (c))!= nullptr ) {
		struct dirent *ent;
		while ( (ent = readdir (dir)) != nullptr ) {
			//TR<<"Loading PSSM from file "<<target_path+ent->d_name<<std::endl;
			utility::io::izstream data( target_path+ent->d_name);
			if ( !data ) {
				continue;
			}
			/*get segment name from file name*/
			string segment_name(target_path+ent->d_name);
			// Remove directory if present.
			// Do this before extension removal incase directory has a period character.
			const size_t last_slash_idx = segment_name.find_last_of("\\/");
			if ( std::string::npos != last_slash_idx ) {
				segment_name.erase(0, last_slash_idx + 1);
			}

			// Remove extension if present.
			const size_t period_idx = segment_name.rfind('.');
			if ( std::string::npos != period_idx ) {
				segment_name.erase(period_idx);
			} else {
				continue;
			}
			// TR<<"segment name:"<<segment_name<<std::endl;
			if ( (segment_name.compare("") == 0)||(segment_name.compare(".") == 0) ) { //hack to avoid trying to read . or.. when looking into directory
				continue;
			}
			read_profile(target_path+ent->d_name,segment_name);
		}
		closedir (dir);
	} else {
		/* could not open directory */
		utility_exit_with_message( "Directory path" + target_path+"was not found" );

	}
}

void //this version of read_profile is an older one and it allows the user to enter it's own PSSM files, instead of using the ones in Rosetta_DB
SpliceSegment::read_profile( string const & file_name, string const & segment_name ){
	SequenceProfileOP new_seqprof( new SequenceProfile );
	utility::file::FileName const fname( file_name );
	//TR<<"Segment name: "<<segment_name<<",reading sequence profile from "<<file_name<<std::endl;
	utility::io::izstream data( fname );
	if ( !data ) {
		utility_exit_with_message( "PSSM file not found " + file_name );
	}
	new_seqprof->read_from_file( fname );
	//TR<< "SpliceSegmentsSeqProf:"<< new_seqprof->prof_row(1)<<std::endl;
	//TR<< "SpliceSegmentspprobabilty:"<< new_seqprof->probability_row(1)<<std::endl;
	sequence_profile_.insert( pair< string, SequenceProfileOP >( segment_name, new_seqprof ) );
}

SequenceProfileOP
SpliceSegment::get_profile( std::string const & segment_name ){

	return sequence_profile_[ segment_name ];
}

void
SpliceSegment::add_pdb_profile_pair( std::string const & pdb, string const & profile_name ){
	pdb_to_profile_map_.insert( pair< string, string >( pdb, profile_name ) );
}

/// @brief this is the most useful method, returning the sequence profile according to a pdb file name
SequenceProfileOP
SpliceSegment::pdb_profile( std::string const & pdb_name ){
	//TR<<"size of sequence_Profile_ is:"<<sequence_profile_.size()<<std::endl;
	//TR<<"pdb to profile map is:"<<pdb_name<<":"<<pdb_to_profile_map_[pdb_name]<<std::endl;
	return sequence_profile_[ pdb_to_profile_map_[ pdb_name ] ];
}

core::sequence::SequenceProfileOP
SpliceSegment::sequence_profile( string const & profile_name ){
	return sequence_profile_[ profile_name ];
}

/// @brief generate one long sequence profile out of several sequence-profile segments.
SequenceProfileOP
concatenate_profiles( utility::vector1< SequenceProfileOP > const & profiles, utility::vector1< std::string > segment_names_ordered, std::string H3_seq){
	SequenceProfileOP concatenated_profile( new SequenceProfile );
	core::Size current_profile_size( 0 );
	core::Size current_segment_name( 1 );
	//TR<<"H3 seq "<<H3_seq<<std::endl;
	SequenceProfileOP H3_profile( new SequenceProfile );
	//TR<<"sequence of H3:"<<H3_seq<<std::endl;
	if ( H3_seq!="" ) {
		core::sequence::Sequence H3seq;
		H3seq.sequence(H3_seq);
		H3_profile->generate_from_sequence(H3seq);
		//TR<<"Size of H3 profile is: "<<H3_profile->size()<<std::endl;
	}


	for ( SequenceProfileCOP prof : profiles ) {
		TR<<"now adding profile of segment "<< segment_names_ordered[ current_segment_name]<<std::endl;
		//find H3 seq and constract a new PSSM from given seqeunce
		bool first_pass=true;
		for ( core::Size pos = 1; pos <= prof->size(); ++pos ) {
			current_profile_size++;
			if ( segment_names_ordered[ current_segment_name]=="H3" && pos>3 && first_pass ) {
				for ( core::Size H3_pos = 1; H3_pos <= H3_profile->size(); ++H3_pos ) {
					concatenated_profile->prof_row( H3_profile->prof_row( H3_pos ), current_profile_size );
					TR<<"Res: "<<current_profile_size<<",prof:"<<H3_profile->prof_row(H3_pos)<<std::endl;
					if ( basic::options::option[ basic::options::OptionKeys::out::file::use_occurrence_data ].value() ) {
						//this option is by default false. If set to true then we read the propensity matrix alongside the pssm matrix. T
						utility::vector1< core::Real > occurence_row (20,0); //H3 lopps apexes do not have occurence data, so I just insert a vector of 0's, gideonla 16/08/14
						concatenated_profile->probabilty_row(occurence_row, current_profile_size );
					}
					current_profile_size++;
				}
				first_pass=false;
			}
			TR<<"Res: "<<current_profile_size<<",prof:"<<prof->prof_row(pos)<<std::endl;
			concatenated_profile->prof_row( prof->prof_row( pos ), current_profile_size );
			concatenated_profile->alphabet(prof->alphabet()); //Save the order of the AA's in the prof_row
			if ( basic::options::option[ basic::options::OptionKeys::out::file::use_occurrence_data ].value() ) {
				//this option is by default false. If set to true then we read the propensity matrix alongside the pssm matrix. T
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
SpliceSegment::read_pdb_profile( std::string const & file_name ){
	utility::io::izstream data( file_name );
	if ( !data ) {
		utility_exit_with_message( "File not found " + file_name );
	}
	//TR<<"Loading pdb profile pairs from file "<<file_name<<std::endl;
	string line;
	while ( getline( data, line ) ) {
		istringstream line_stream( line );
		while ( !line_stream.eof() ) {
			std::string pdb, profile;
			line_stream >> pdb >> profile;
			pdb_to_profile_map_.insert( pair< string, string >( pdb, profile ) );
			//TR<<"Loading pdb-profile pair: "<<pdb<<" "<<profile<<std::endl;
		}
	}
}

void
SpliceSegment::add_profile( string const & pdb_name,  SequenceProfileOP const prof){
	sequence_profile_.insert( pair< string, SequenceProfileOP >( pdb_name, prof ) );
}

void
SpliceSegment::read_pdb_profile_file( string const & Protein_family_path, string const segment){

	//This is to get the path to the database
	utility::io::izstream data;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Get all PDB_profile_Match files from raltive path
	const std::string target_path = basic::database::full_name( Protein_family_path+"pdb_profile_match/" );


	DIR *dir;
	const char * c =target_path.c_str();
	if ( (dir = opendir (c))!= nullptr ) {
		string fileName(target_path+"pdb_profile_match."+segment);
		utility::io::izstream data( fileName );
		if ( !data ) {
			utility_exit_with_message( "File not found " + fileName );
		}
		std::string line;
		// commenting out with gideon's permission -- producing silly changes in integration tests due to path differences.
		//TR<<"Loading pdb profile pairs from file "<<fileName<<std::endl;
		while ( getline( data, line ) ) {
			istringstream line_stream( line );
			while ( !line_stream.eof() ) {
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

/// @brief Reads from given file the H3 sequences from all PDBs in the db
std::map< std::string, std::string>
read_H3_seq( std::string const & Protein_family_path){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::map< std::string, std::string>pdb_to_H3_seq_map_;
	pdb_to_H3_seq_map_.clear();
	const std::string H3_seq_file = basic::database::full_name( Protein_family_path+"pssm/H3/H3_seq" );
	utility::io::izstream data( H3_seq_file );
	if ( !data ) {
		TR<<"!!Can't open H3 seq file: "<<H3_seq_file<<". This is probably an error. Using PSSM files instead!!"<<std::endl;
		return pdb_to_H3_seq_map_;//if file not found then return empty map.
	}
	std::string line;
	while ( getline( data, line ) ) {
		istringstream line_stream( line );
		while ( !line_stream.eof() ) {
			std::string pdb, seq;
			line_stream >> pdb >> seq;
			pdb_to_H3_seq_map_.insert( pair< string, string >( pdb, seq ) );
			//TR<<"Loading pdb-seq pair: "<<pdb<<" "<<seq<<std::endl;
		}
	}
	return pdb_to_H3_seq_map_;
}
///@brief in cases where we have only one sequence per pdb we store the sequences in a seperate map
///we can later use them to create pssms using blosum
std::map< std::string, std::string>
read_pdb_seq( std::string const & Protein_family_path){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::map< std::string, std::string> pdb_seq_map;
	pdb_seq_map.clear();
	const std::string seq_file = basic::database::full_name( Protein_family_path );
	utility::io::izstream data( seq_file );
	if ( !data ) {
		TR<<"!!Can't open seq file: "<<seq_file<<". This is probably an error. Using PSSM files instead!!"<<std::endl;
		return pdb_seq_map;//if file not found then return empty map.
	}
	std::string line;
	while ( getline( data, line ) ) {
		istringstream line_stream( line );
		while ( !line_stream.eof() ) {
			std::string pdb, seq;
			line_stream >> pdb >> seq;
			pdb_seq_map.insert( pair< string, string >( pdb, seq ) );
			//TR<<"Loading pdb-seq pair: "<<pdb<<" "<<seq<<std::endl;
		}
	}
	return pdb_seq_map;
}


} //splice
} //protocols
