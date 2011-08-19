// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <protocols/ScoreMap.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>

// Auto-header: duplicate removed #include <core/io/pdb/pose_io.hh>

#include <numeric/random/random.hh>

#include <protocols/ddg/ddGData.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>


// AUTO-REMOVED #include <sstream>

// AUTO-REMOVED #include <basic/Tracer.hh>

numeric::random::RandomGenerator RG(15431); // <- Magic number, do not change it!!!


namespace protocols {
namespace ddG{

ddGData::ddGData(){
	file_to_read_ = "";
	curr_mut_filename = "";
	curr_wt_filename = "";
	experimental = 0.0;
	curr_mut_read = true;
	curr_wt_read = true;
}

ddGData::ddGData(std::string filename){
	file_to_read_ = filename;
	curr_mut_filename = "";
	curr_wt_filename = "";
	experimental = 0.0;
	curr_mut_read = true;
	curr_wt_read = true;
}

//returns false if file hasn't been opened
bool ddGData::end(){
	if(!inputstream.is_open() && (file_to_read_.compare("") != 0)){
		inputstream.open(file_to_read_.c_str());
	}
	return (inputstream.eof());
}

void ddGData::get_next_filenames(){
	if(!inputstream.is_open()){
		inputstream.open(file_to_read_.c_str());
	}
	if(inputstream.is_open() && !inputstream.eof() && curr_mut_read && curr_wt_read){
		inputstream >> curr_wt_filename >> curr_mut_filename >> experimental;
	}
	curr_mut_read = false;
	curr_wt_read = false;
}

utility::vector1<pose::Pose> ddGData::read_mut_data()
{
	utility::vector1<pose::Pose> curr_mut;
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	if(!inputstream.eof() && !curr_mut_read){
		curr_mut_read = true;
	}
	utility::file::FileName mut(curr_mut_filename);
	if(mut.ext().compare("list") == 0){
		std::ifstream istream;
		istream.open(curr_mut_filename.c_str());
		if(istream.is_open()){
			while(!istream.eof()){
				pose::Pose newpose;
				std::string curr_filename;
				istream >> curr_filename;
				if(curr_filename.compare("") != 0){
					std::cout << "opening file: " << curr_filename << std::endl;
					core::import_pose::pose_from_pdb(newpose,(*rsd_set),curr_filename,false);
					curr_mut.push_back(newpose);
				}
			}
		}
	}else if(mut.ext().compare("out") == 0){
		core::io::silent::SilentFileData sfd;
		sfd.set_filename(curr_mut_filename);
		utility::vector1<std::string> tags = sfd.tags();
		for(unsigned int i=1;i <= tags.size(); i++){
			core::io::silent::SilentStructOP ss = sfd[tags[i]];
			pose::Pose newpose;
			ss->fill_pose(newpose,(*rsd_set));
			curr_mut.push_back(newpose);
		}
	}else{
		//extension not recognized
	}
	return curr_mut;
}

utility::vector1<pose::Pose> ddGData::read_wt_data(){
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	utility::vector1<pose::Pose> curr_wt;
	if(!inputstream.eof() && !curr_wt_read){
		curr_wt_read = true;
	}
	utility::file::FileName wt(curr_wt_filename);
	if(wt.ext().compare("list") == 0){
		std::ifstream istream;
		istream.open(curr_wt_filename.c_str());
		if(istream.is_open()){
			while(!istream.eof()){
				pose::Pose newpose;
				std::string curr_filename;
				istream >> curr_filename;
				if(curr_filename.compare("") != 0){
					core::import_pose::pose_from_pdb(newpose,(*rsd_set),curr_filename,false);
					curr_wt.push_back(newpose);
				}
			}
		}
	}else if(wt.ext().compare("out") == 0){
		core::io::silent::SilentFileData sfd;
		sfd.set_filename(curr_wt_filename);
		utility::vector1<std::string> tags = sfd.tags();
		for(unsigned int i=1;i <= tags.size(); i++){
			core::io::silent::SilentStructOP ss = sfd[tags[i]];
			pose::Pose newpose;
			ss->fill_pose(newpose,(*rsd_set));
			curr_wt.push_back(newpose);
		}
	}else{
		//extension not recognized
	}
	return curr_wt;
}

core::Real ddGData::read_exp_data(){
	return experimental;
}


PositionDdGInfo::PositionDdGInfo(
	core::Size seqpos,
  core::chemical::AA wt_aa
)
	: seqpos_(seqpos),
		wt_aa_(wt_aa)
{}

PositionDdGInfo::~PositionDdGInfo(){}

/// @details doesn't check whether something for this mutation already exists,
/// so previously added stuff will be overwritten
void
PositionDdGInfo::add_mutation_ddG(
   core::chemical::AA aa,
   core::Real ddG
)
{
	std::map< core::chemical::AA, core::Real >::iterator map_it = mutation_ddGs_.find( aa );
	if( map_it == mutation_ddGs_.end() ){
		mutation_ddGs_.insert( std::pair< core::chemical::AA, core::Real >(aa, ddG) );
	}
	else map_it->second = ddG;
}

const std::map< core::Size, PositionDdGInfoOP >
read_ddg_predictions_file( std::string filename )
{

	utility::io::izstream data( filename.c_str() );
	if ( !data ) utility_exit_with_message("File "+filename+"could not be opened.");
	std::map< core::Size, PositionDdGInfoOP > to_return;

	while( !data.eof() ) {

		std::string line;
		getline(data,line);
		utility::vector1< std::string > tokens; tokens.push_back("");
		tokens = utility::split(line);
		if( tokens.size() < 2 ) continue; //ddg predictions output file has empty lines :(
		std::string mutstring( tokens[2] );
		core::Size mutaaloc( mutstring.length() - 1 );
		std::string wt_res = mutstring.substr(0,1);
		std::string mut_res = mutstring.substr(mutaaloc, 1 );
		std::string string_seqpos( mutstring.substr(1, mutaaloc - 1 ) );
		core::Size mut_seqpos( atoi( string_seqpos.c_str() ));
		core::Real ddG( atof(tokens[3].c_str() ) );

		//std::cerr << "Reading ddg out file, at pos " << mut_seqpos << " wt " << wt_res << " turned to " << mut_res << " with ddG of " << ddG << std::endl;

		core::chemical::AA wt_aa( core::chemical::aa_from_oneletter_code( wt_res[0] ) );
		core::chemical::AA mut_aa( core::chemical::aa_from_oneletter_code( mut_res[0] ) );

		std::map< core::Size, PositionDdGInfoOP >::iterator muts_it( to_return.find( mut_seqpos ));
		if( muts_it == to_return.end() ){
			PositionDdGInfoOP pos_info = new PositionDdGInfo( mut_seqpos, wt_aa );
			to_return.insert( std::pair< core::Size, PositionDdGInfoOP >( mut_seqpos, pos_info ) );
			muts_it = to_return.find( mut_seqpos );
		}
		muts_it->second->add_mutation_ddG( mut_aa, ddG );
	}
	data.close();
	return to_return;
}

} //namespace ddG
} //namespace protocols
