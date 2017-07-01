// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/forge/remodel/RemodelData.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Rosetta Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/file/FileName.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <core/chemical/ChemicalManager.hh> // for ncaa handling
#include <core/chemical/ResidueTypeSet.hh> // for ncaa handling
#include <core/chemical/ResidueType.hh> // for ncaa handling

#include <protocols/forge/remodel/RemodelData.hh>
//#include <protocols/sparta/PDB.hh>
#include <core/pose/PDBInfo.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
// C++ Headers
#include <vector>
#include <string>


static THREAD_LOCAL basic::Tracer TR_REMODEL( "protocols.forge.remodel.RemodelData" );

using namespace basic::options::OptionKeys;
using namespace basic::options;
using namespace core;

namespace protocols {
namespace forge {
namespace remodel {

///
/// @brief
/// constructor
///
RemodelData::RemodelData()
{
	has_design_info_ = false;
	design_neighbor = false;
	auto_design = false;
	natro_movemap_.set_chi( true );
}

///
/// @brief
/// Splits a string on the passed in delimeter and places the tokens in the passed in vector.
/// Isn't there a function in utility that can do this already?
///
void RemodelData::splitString( std::string str, std::string delim, std::vector< std::string > & results ) {
	int cutAt;
	while ( (cutAt = (int)str.find_first_of(delim)) != int(str.npos ) ) {
		if ( cutAt > 0 ) {
			results.push_back(str.substr(0,cutAt));
		}
		str = str.substr(cutAt+1);
	}
	if ( str.length() > 0 ) {
		results.push_back( str );
	}
}

/// @brief
/// Parses the blueprint text
void RemodelData::getLoopsToBuildFromBlueprint( std::string text_blueprint ) {
	using namespace ObjexxFCL;
	if ( text_blueprint == "" ) {
		TR_REMODEL << "No blueprint data given!" << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}

	std::stringstream data(text_blueprint);
	std::string line;

	// extension management
	//std::string ext_ss_buffer;

	std::ostringstream oss; // for resfile parsing
	//std::ostringstream oss_switch; // for detecting design info
	oss << "NATRO" << std::endl; // preseve non designed to starting rotamer
	oss << "start" << std::endl; // mark start for resfile parser

	// first getting the total line in the file, needed for repeat Resfile processing
	int index = 0;
	while ( getline( data, line ) ) {
		index++;
	}
	int length = index;

	// reset file and counter
	data.clear();
	//data.seek_beg();
	data.seekg(0,std::ios::beg);
	index = 1;

	bool mark_start = false;
	while ( getline( data, line ) ) {

		std::istringstream line_stream( line );
		std::vector< std::string > split_info;
		std::ostringstream oss_switch; // for detecting design info
		// split line on spaces and store resulting strings in split_info
		splitString( line_stream.str(), " ", split_info );

		// skip comment lines
		if ( split_info[0].at(0) == '#' ) {
			continue;
		}

		// remodel::LineObject is a struct which has several field incl. index, original_index, resname, sstype,
		// design_type, isDesignable, has_constraints, constraint_definition, and aminoAcidList
		protocols::forge::remodel::LineObject line;
		line.index = index;
		line.isDesignable = false; // initialize design by default to false
		line.has_ncaa = false;

		// debug
		//line_stream >> line.original_index >> line.resname >> line.sstype >> skip;
		//this->blueprint.push_back(line);
		//
		index++;

		// could have initialized blueprint after the split, oh well...
		std::istringstream( split_info[0] ) >> line.original_index;
		if ( line.original_index != 0 && mark_start == false ) {
			pdb_start = line.original_index;
			mark_start = true;
		}
		if ( line.original_index != 0 ) {
			pdb_stop = line.original_index;
		}

		line.resname = split_info[1];
		line.sstype = split_info[2]; // could be 'H', 'L', 'E', 'D', or 'I", which stand for helix, loop, extended, random or inserted
		//Feb 2016. Extended to handle ABEGO and also D-amino acids '1', '2', '3'.

		// do some input file syntax checking
		// error checking, disallow '#' in columns 1 and 2
		if ( split_info[1].at(0) == '#' || split_info[2].at(0) == '#' ) {
			std::ostringstream err_message;
			err_message << "ERROR: comment marker '#' cannot be in residue or ss column at line:\n";
			err_message << line_stream.str();
			utility::exit( __FILE__, __LINE__, err_message.str() );
		}

		// if the blueprint line has something like "20 F . ALLAA" then the no. of tokens will be > 3 and ALLAA is the design info
		// need to also figure out whether these lines need to be part of the resfile
		if ( split_info.size() > 3 ) { // has design info

			// skip comments at end of line
			if ( split_info[3].at(0) == '#' ) {
				continue;
			}

			// check for the presence of special tokens, like CST, DM_, DS_, in the blueprint file
			for ( std::vector< std::string >::iterator it = split_info.begin(), end = split_info.end(); it != end; ++it ) {

				TR_REMODEL << "Found token " << *it << " on blueprint line '" << line.index << "'" << std::endl;

				if ( (*it).substr(0,3) == "CST" ) {
					TR_REMODEL << "constraint found " << *it <<  std::endl;
					line.has_constraints=true;
					line.constraint_definition.push_back(*it);
				}
				if ( (*it).substr(0,3) == "DM_" ) {
					disulfMobileRange.push_back(line.index);
					TR_REMODEL << "Added line " << line.index << " to disulfMobileRange while reading blueprint file." << std::endl;
					if ( disulfMobileRange.size() > 2 ) {
						std::ostringstream err_message;
						err_message << "ERROR: Disulfide mobile range assigment contains " <<  disulfMobileRange.size() << " elements." << std::endl;
						utility::exit( __FILE__, __LINE__, err_message.str() );
					}
				}
				if ( (*it).substr(0,3) == "DS_" ) {
					disulfLandingRange.push_back(line.index);
					if ( disulfLandingRange.size() > 2 ) {
						std::ostringstream err_message;
						err_message << "ERROR: Disulfide landing range assigment contains " <<  disulfLandingRange.size() << " elements." << std::endl;
						utility::exit( __FILE__, __LINE__, err_message.str() );
					}
				}

			}

			std::map< std::string, pack::task::ResfileCommandOP > resfile_command_map = pack::task::create_command_map();

			// need to figure out if the line is to be part of Resfile
			bool design_info = false;
			int info_count = (int)split_info.size() - 3; // only keep count past 3
			//std::cout << "split info count: " << (int) split_info.size()  << std::endl;

			for ( int i = 3; i< (int)split_info.size(); i++ ) {
				if ( split_info[i].substr(0,3) == "CST" || split_info[i].substr(0,3) == "DM_" || split_info[i].substr(0,3) == "DS_" ) {
					//std::cout << " split info: " << split_info[i].substr(0,3) << std::endl;
					info_count--;
				}
			}
			if ( info_count > 0 ) {
				design_info = true;
			}


			if ( design_info ) {
				if ( option[OptionKeys::run::chain].user() ) {
					std::string const chain(option[OptionKeys::run::chain] );
					oss << line.index << " " << chain << " " ;
					TR_REMODEL << "Do remodeling on chain " << chain << std::endl;
				} else {

					/* JAB - can't rely on using a single PDB file. It might be CIF/etc, so for now, we can't use this logic!

					// For the output resfile, 'chain' is defined as the the first chain in the input pdb!
					// If a list is used, we assume that all files in that list have the same amount and lengthes of chains.
					// This is because otherwise the blueprint file wouldn't be universal for all those files.
					// So, for now we take the chain ID of the first chain in the first structur ein the list
					// If this turns out to be not flexible enough, one would need to get the current pose from the
					// job distributer. [Sebastian Rämisch and Jared Adolf-Bryfogle, Jan 2016.]

					// Take the first chain in the input pdb file
					core::pose::Pose inputPose;
					std::string input_pdb;
					if ( option[OptionKeys::in::file::s].user()  )
					input_pdb = option[ in::file::s ]()[1];
					if ( option[OptionKeys::in::file::l].user()  ) {
					utility::vector1<utility::file::FileName> list(basic::options::option[ in::file::l ]() );
					utility::vector1<std::string> files;
					for ( unsigned int h=1; h<=list.size(); h++ ) {
					utility::io::izstream pdbs(list[h]);
					std::string fname;
					while ( pdbs >> fname ) {
					files.push_back(fname);
					}
					}
					input_pdb = files[1];
					}
					import_pose::pose_from_pdb( inputPose, input_pdb );
					char chain(inputPose.pdb_info()->chain(1));
					oss << line.index << " " << chain << " " ;
					TR_REMODEL << "-run::chain not found. No chain specified." << std::endl;
					TR_REMODEL << "Use first chain in input pdb: " << std::endl;
					TR_REMODEL << "... do remodeling on chain " << chain << std::endl;
					*/



					oss << line.index << " A " ; //default to chain A, Jan 2016

				}
			}

			bool pickaa = false;
			for ( int i = 3; i< (int)split_info.size();  i++ ) {

				if ( split_info[i].substr(0,3) != "CST" && split_info[i].substr(0,3) != "DM_" && split_info[i].substr(0,3) != "DS_" && split_info[i].substr(0,3) != "EMP" ) {
					oss << split_info[i] << " " ;
					oss_switch << split_info[i];
				}
				if ( split_info[i].substr(0,5) == "PIKAA" ) {
					// toggle on manual residue selection switch
					pickaa = true;
					continue;
				}
				if ( split_info[i].substr(0,5) == "EMPTY" ) {
					line.has_ncaa = true;
					// get the list of ncaa
					oss << split_info[i];
					for ( int j = i+1; j < (int)split_info.size() - 1; j = j+2 ) {
						if ( split_info[j] == "NC" ) {
							core::chemical::ResidueTypeSetCOP res_type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
							std::string ncaa_fullname = res_type_set->name_map(split_info[j+1]).name3();
							std::cout << "debugDATA " << ncaa_fullname << std::endl;
							line.ncaaList.push_back(split_info[j+1]);
							oss << " " << split_info[j] << " "  << ncaa_fullname;
							i+=2;
						}
					}
				}
				if ( pickaa ) { // the column following PIKAA
					for ( int j=0; j < (int)split_info[i].size(); ++j ) { // only find string element size

						//char one_letter_name(core::chemical::aa_from_oneletter_code(aa));
						char one_letter_name = split_info[i].substr(j,1).c_str()[0];
						chemical::AA aa = core::chemical::aa_from_oneletter_code( one_letter_name );
						TR_REMODEL << "  design position " << line.index << " to " << one_letter_name << " " << aa << std::endl;

						line.aminoAcidList.push_back(aa);
					}
					pickaa = false; // turns it right off so doesn't get other columns
				}

				if ( split_info[i] == "NATRO" ) {
					TR_REMODEL << "NATRO movemap setup: turning off chi move for refinement stage: " << line.index << std::endl;
					natro_movemap_.set_chi(line.index, false);
				}
			}

			// add a newline to the stream if there was design information on this line of the blueprint file
			if ( design_info ) {
				oss << std::endl;
			}

			// process repeats, pretty dangerous, as this only hacks the resfile string
			// but not making duplicates in the blueprint held by RemodelData
			if ( option[OptionKeys::remodel::repeat_structure].user() ) {
				for ( Size rep = 1; rep <(Size)option[OptionKeys::remodel::repeat_structure]; rep++ ) {
					//chain defined by option, no chain by default
					if ( option[OptionKeys::run::chain].user() ) {
						std::string const chain(option[OptionKeys::run::chain]);
						oss << line.index + length*rep << " " << chain << " " ;
					} else {
						oss << line.index + length*rep << " _ " ;
					}
					for ( Size i = 3; i< split_info.size(); i++ ) {
						if ( split_info[i].substr(0,3) != "CST" ) {
							oss << split_info[i] << " " ;
						}
					}
					oss << std::endl;
				}
			}

			// find out that there's info other than CST and turn on manual modes
			if ( oss_switch.str() != "" ) {
				//TR_REMODEL << "oss_switch: " << oss_switch.str() << std::endl;
				has_design_info_ = true;
			}


			parsed_string_for_resfile = oss.str();

			//TR_REMODEL << "manual design overwrite position: " << line.index << std::endl;
			//this->design_mode = 3; //default manual mode
			/*if (option[OptionKeys::Design::design_neighbors] ) {
			// fully manual design mode automatically switched on when you assign residues by hand
			this->design_mode = 4;
			}
			if(option[OptionKeys::Design::neighbor_repack] ) {
			// bc repack neigbors
			this->design_mode = 5;
			}
			*/
			line.isDesignable = true;
			// check to make sure the design type (i.e. PIKAA, NATRO, NATAA) is one of the recognized packer tokens
			line.design_type = split_info[3];
			if ( !resfile_command_map[ line.design_type ] ) {
				TR_REMODEL.Warning << "unknown packer token: " << line.design_type << std::endl;
			}

			//debug
			//TR_REMODEL << resfile_command_map[split_info[3]] << " resfile command map to " << split_info[3] << std::endl;

		}

		blueprint.push_back(line);

	}

	if ( parsed_string_for_resfile != "" ) {
		TR_REMODEL << "Resfile to be used for design:\n" << parsed_string_for_resfile << std::endl;
	}

	if ( has_design_info_ ) {
		design_mode = 3; // what is design mode 3?!?
	}

	// process blueprint to initialize all the needed strings/vectors
	// iterate over all the LineObject objects and save the sequence and ss.
	bool hle_abego_mode = false;
	std::vector< protocols::forge::remodel::LineObject >::iterator iter, end;
	for ( iter = this->blueprint.begin(), end = this->blueprint.end(); iter != end; ++iter ) {
		// sequence and ss are class member variables
		sequence.append( iter->resname );
		if ( iter->sstype.size()==1 ) {
			//Feb 2016.  handling LD types
			if ( iter->sstype.compare( "1" )  == 0 ) { // 0 means match
				ss.append( "H" );
				LD_types.append( "D" );
			} else if ( iter->sstype.compare( "2" ) == 0 ) {
				ss.append( "E" );
				LD_types.append( "D" );
			} else {
				ss.append( iter->sstype ); //std case
				LD_types.append( "L" ); //std L-amino acid case
			}
		} else if ( iter->sstype.size()==2 ) {
			hle_abego_mode=true;
			char tmp_ss = iter->sstype[0];
			char tmp_abego = iter->sstype[1];
			ss.append(1,tmp_ss);
			abego.append(1,tmp_abego);
			LD_types.append( "L" );
			if ( !(tmp_ss=='H'||tmp_ss=='L'||tmp_ss=='E'||tmp_ss=='D') ) {
				utility_exit_with_message("First SS-term is:" + string_of(tmp_ss)+ " but must be either H,L,E or D if you want it ignored");
			}
			if ( !(tmp_abego=='A'||tmp_abego=='B'||tmp_abego=='E'||tmp_abego=='G'||tmp_abego=='O'||tmp_abego=='D') ) {
				utility_exit_with_message("Second SS-term is:" +string_of(tmp_abego)+" but must be either A,B,E,G,O or D if you want it ignored");
			}
		}
		if ( iter->sstype.size()==1 && hle_abego_mode==true ) {
			utility_exit_with_message("Blueprint style must be either all abego,HLE or HLE followed by ABEGO but not both");
		}
		if ( iter->sstype.size()==3 ) {
			utility_exit_with_message("Blueprint must have no more than 2 ssTYPES. Ex. HA");
		}
	}
	if ( hle_abego_mode==false ) {
		translateDSSP_ABEGO(this->ss, this->abego);
	}

	TR_REMODEL << " sequence (based on blueprint): " << std::endl << " " << sequence << std::endl;
	TR_REMODEL << " ss (based on blueprint): " << std::endl << " " << ss << std::endl;
	TR_REMODEL << " ABEGOtype: " << std::endl << abego << std::endl;
	TR_REMODEL << " LD_types: " << std::endl << LD_types << std::endl;

}


/// @brief Reads in the blueprint file and passes the text data to the blueprint file parser
void RemodelData::getLoopsToBuildFromFile( std::string filename) {
	using namespace ObjexxFCL;
	if ( filename == "" ) {
		TR_REMODEL << "No blueprint file given!" << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__ );
	}

	std::stringstream data;
	std::string line;

	utility::io::izstream file_data( filename.c_str() );
	while ( getline( file_data, line ) ) {
		data << line << std::endl;
	}

	if ( !data ) {
		TR_REMODEL << "Can't open blueprint file " << filename << std::endl;
		utility::exit(EXIT_FAILURE, __FILE__, __LINE__);
	}

	getLoopsToBuildFromBlueprint(data.str());
}

///
/// @brief
/// Updates the dssp_updated_ss vector with secondary structure information. Uses the information obtained from
/// DSSP unless the secondary structure was specified in the blueprint file.
///
void RemodelData::updateWithDsspAssignment( ObjexxFCL::FArray1D_char & dsspSS ) {

	for ( int i = 0; i < (int)ss.size(); i++ ) {
		int idx = blueprint[i].original_index;
		char const * ss_chars = ss.c_str();
		if ( ss_chars[i] != '.' ) {
			dssp_updated_ss.append( 1, ss_chars[i] );
		} else {
			dssp_updated_ss.append( 1, dsspSS(idx) );
		}
	}
	// turn upper case if not already so
	std::transform( dssp_updated_ss.begin(), dssp_updated_ss.end(), dssp_updated_ss.begin(), ::toupper );
	TR_REMODEL << " dssp_updated_ss: lengths = "  << dssp_updated_ss.length() << std::endl << dssp_updated_ss << std::endl;
}


void RemodelData::translateDSSP_ABEGO( std::string & ss, std::string & abego ) {

	size_t found_idx;
	bool abego_switch = false;
	found_idx = ss.find_first_of("abgoABGO"); // E is shared so only ABGO for mapping
	if ( found_idx == std::string::npos ) {
		TR_REMODEL << "SS based assignment found" << std::endl; // in case of only E assignment, treat it as DSSP
	} else if ( found_idx != std::string::npos ) {
		abego_switch = true;
		TR_REMODEL << "ABEGO based assignment found" << std::endl;
	}
	std::string trans_ss;

	if ( abego_switch ) { //need to make a new string with DSSP assignment and swap
		//found_idx = ss.find_first_of("abegoABEGO"); // this substitution use all 5 regions
		for ( core::Size idx = 0; idx < ss.length(); idx++ ) {
			if ( ss[idx] == 'A' || ss[idx] == 'a' ) {
				trans_ss.push_back('D');
			} else if ( ss[idx] == 'B' || ss[idx] == 'b' ) {
				trans_ss.push_back('D');
			} else if ( ss[idx] == 'E' || ss[idx] == 'e' ) {
				trans_ss.push_back('D');
			} else if ( ss[idx] == 'G' || ss[idx] == 'g' ) {
				trans_ss.push_back('D');
			} else if ( ss[idx] == 'O' || ss[idx] == 'o' ) {
				trans_ss.push_back('D');
			} else if ( ss[idx] == '.' ) {
				trans_ss.push_back('.');
				ss[idx]= 'X';
			} else { // could have other characters like I, or D so leave them alone
				trans_ss.push_back(ss[idx]);
			}
		}
		trans_ss.swap(ss);
	}
	abego = trans_ss;

}

///
/// @brief
/// If users are trying to do domain insertion with remodel (by specifying the insert_segment_from_pdb command-line
/// option), then this function is called to read in the PDB file for that segment.
///
void RemodelData::collectInsertionPose(){

	using namespace core;

	import_pose::pose_from_file( insertPose,option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb] , core::import_pose::PDB_file);
	insertionSize = (int)insertPose.size();
	TR_REMODEL << "insertionSize: " << insertionSize << std::endl;

	scoring::dssp::Dssp dssp( insertPose );
	ObjexxFCL::FArray1D_char dsspSS( (int)insertPose.size() );
	dssp.dssp_reduced( dsspSS );
	TR_REMODEL << "dsspSS.size(): " << dsspSS.size() << std::endl;
	for ( int i = 1; i <= (int)dsspSS.size(); i++ ) {
		insertionSS.push_back( dsspSS(i) );
	}
	TR_REMODEL << "insertion SS: " << insertionSS << " " << std::endl;
	//hack for jack
	/*
	import_pose::pose_from_pdb( insertPose2,option[OptionKeys::remodel::domainFusion::insert_segment2_from_pdb] );
	insertion2Size = (int)insertPose2.size();
	TR_REMODEL << "insertionSize: " << insertion2Size << std::endl;

	scoring::dssp::Dssp dssp2( insertPose2 );
	ObjexxFCL::FArray1D_char dssp2SS( (int)insertPose2.size() );
	dssp2.dssp_reduced( dssp2SS );
	TR_REMODEL << "dsspSS.size(): " << dssp2SS.size() << std::endl;
	for ( int i = 1; i <= (int)dssp2SS.size(); i++ ) {
	insertion2SS.push_back( dssp2SS(i) );
	}
	TR_REMODEL << "insertion2 SS: " << insertion2SS << " " << std::endl;
	*/


}

} //namespace remodel
} //namespace forge
} //namespace protocols

