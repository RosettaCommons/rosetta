// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyNumberingParser.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/AntibodyNumberingParser.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>

#include <utility>
#include <utility/string_util.hh>
#include <utility/py/PyAssert.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include <sstream>


namespace protocols {
namespace antibody {

static basic::Tracer TR( "protocols.antibody.AntibodyNumberingParser" );

//using namespace boost;
using utility::vector1;

AntibodyNumberingParser::AntibodyNumberingParser(AntibodyEnumManagerCOP enum_manager):
	enum_manager_(std::move(enum_manager)),
	numbering_database_directory_("sampling/antibodies/numbering_schemes"),
	scheme_file_("numbering_scheme_definitions.txt"),
	cdr_definition_file_("cdr_definitions.txt")
{}

AntibodyNumberingParser::~AntibodyNumberingParser()= default;

void
AntibodyNumberingParser::check_path(std::string const & numbering_file_path) const {
	using namespace std;
	ifstream check( numbering_file_path.c_str(), ifstream::in);
	if ( check.good() ) { return;}
	else {
		utility_exit_with_message("Antibody Numbering file does not exist.  Cannot continue. "+numbering_file_path);
	}
}


AntibodyNumbering
AntibodyNumberingParser::get_antibody_numbering(AntibodyNumberingSchemeEnum const numbering_scheme, CDRDefinitionEnum const cdr_definition){


	//Clear anything.
	cdr_definitions_defined_.clear();
	cdr_definitions_defined_using_.clear();
	schemes_defined_.clear();
	std::string cdr_def_file_path = basic::database::full_name(numbering_database_directory_+"/"+cdr_definition_file_);
	std::string scheme_def_file_path = basic::database::full_name(numbering_database_directory_ + "/"+scheme_file_);

	AntibodyNumbering numbering;


	numbering.numbering_scheme = numbering_scheme;
	numbering.cdr_definition = cdr_definition;

	read_numbering_scheme_file(scheme_def_file_path, numbering);
	read_cdr_definition_file(cdr_def_file_path, numbering);
	//debug_print(numbering);
	return numbering;
}

void
AntibodyNumberingParser::read_numbering_scheme_file(std::string const & file_path, AntibodyNumbering& numbering) {
	check_path(file_path);
	std::ifstream numbering_file(file_path.c_str());
	PyAssert((numbering_file.is_open()), "Unable to open scheme transform file.");

	std::string line;

	while ( getline(numbering_file, line) ) {

		//Skip any comments + empty lines
		if ( utility::startswith(line, "#") ) {
			continue;
		}
		if ( utility::startswith(line, "\n") ) {
			continue;
		}

		utility::trim(line, "\n"); //Remove trailing line break

		vector1< std::string > lineSP = utility::string_split_multi_delim(line); //Split on space or tab

		if ( (lineSP.size()) == 0 ) continue;

		if ( lineSP[1] == "DEFINES" ) {
			read_scheme_defines_line(lineSP);
		} else {
			read_scheme_numbering_line(lineSP, numbering);
		}
	}

	core::Size landmarks_defined = numbering.numbering_scheme_transform[schemes_defined_[1]].size();
	for ( core::Size i = 2; i <= schemes_defined_.size(); ++i ) {
		if ( landmarks_defined != numbering.numbering_scheme_transform[schemes_defined_[i]].size() ) {
			utility_exit_with_message("One to one correspondence in scheme transforms is needed.  Please make sure all PDB numbers have equivalent positions in each numbering scheme in scheme file.");

		}
	}
	TR << "Antibody numbering scheme definitions read successfully" << std::endl;
	numbering_file.close();

}

void
AntibodyNumberingParser::read_scheme_defines_line(vector1< std::string > const & lineSP) {

	for ( core::Size i = 2; i <= lineSP.size(); ++i ) {
		if ( enum_manager_->numbering_scheme_is_present(lineSP[i]) ) {
			schemes_defined_.push_back(enum_manager_->numbering_scheme_string_to_enum(lineSP[i]));
		} else {
			utility_exit_with_message("Numbering scheme unrecognized while parsing scheme file: "+lineSP[i]);
		}
	}
}

void
AntibodyNumberingParser::read_scheme_numbering_line(vector1< std::string > const & lineSP, AntibodyNumbering & numbering) const {

	if ( schemes_defined_.size() != lineSP.size() ) {
		utility_exit_with_message("Transform schemes from scheme file do match number of defined schemes.  ");
	}

	for ( core::Size i = 1; i <= lineSP.size(); ++i ) {
		vector1< std::string > raw_landmark = utility::string_split(lineSP[i], ':');
		if ( raw_landmark.size() != 3 ) {
			utility_exit_with_message("Error while reading scheme file.  landmark must define: chain:resnum:insertion code.  If no insertion code, please use a period (.)");
		}

		char insertion_code;
		if ( raw_landmark[3] == "~" ) {
			insertion_code = ' ';
		} else {
			insertion_code = raw_landmark[3].at(0);
		}

		char chain = raw_landmark[1].at(0);
		core::Size resnum; std::stringstream(raw_landmark[2]) >> resnum;
		PDBLandmarkOP landmark( new PDBLandmark(chain, resnum, insertion_code, schemes_defined_[i]) );

		numbering.numbering_scheme_transform[schemes_defined_[i]].push_back(landmark);
	}
}

void
AntibodyNumberingParser::debug_print(AntibodyNumbering& numbering)  {

	//for (core::Size i = 1; i <= numbering.numbering_scheme_transform[Chothia_Scheme].size(); ++i){
	// TR << numbering.numbering_scheme_transform[Chothia_Scheme][i]->get_string()<< " "
	// << numbering.numbering_scheme_transform[AHO_Scheme][i]->get_string() << std::endl;
	//}
	TR << "Scheme: "<< enum_manager_->numbering_scheme_enum_to_string(numbering.numbering_scheme)<<std::endl;
	TR << "Definition:" << enum_manager_->cdr_definition_enum_to_string(numbering.cdr_definition)<< std::endl;
	for ( core::Size i = 1; i <= 6; ++i ) {
		//CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		TR<<"North cdr_start: "<<numbering.cdr_definition_transform[North][i][cdr_start]->get_string() << " "
			<< "Chothia cdr_start:: " << numbering.cdr_definition_transform[Chothia][i][cdr_start]->get_string() << std::endl;
	}
}

PDBLandmarkOP
AntibodyNumberingParser::get_equivalent_landmark(AntibodyNumbering & numbering, const AntibodyNumberingSchemeEnum scheme, PDBLandmark & landmark_to_match) const {


	std::map< AntibodyNumberingSchemeEnum, vector1< PDBLandmarkOP > > ::const_iterator iter( numbering.numbering_scheme_transform.find( landmark_to_match.numbering_scheme() ) );
	if ( iter == numbering.numbering_scheme_transform.end() ) {
		utility_exit_with_message("AntibodyNumberingScheme used to define a cdr definition is not defined in numbering_scheme_definitions file");
	}
	vector1< PDBLandmarkOP > landmarks = numbering.numbering_scheme_transform[scheme];

	//TR<< "Need "<< enum_manager_->numbering_scheme_enum_to_string(landmark_to_match.numbering_scheme())<< " in " << enum_manager_->numbering_scheme_enum_to_string(numbering.numbering_scheme) << std::endl;
	for ( core::Size i = 1; i <= landmarks.size(); ++i ) {
		PDBLandmarkOP landmark =landmarks[i];
		if ( *landmark == landmark_to_match ) {

			PDBLandmarkOP new_landmark = numbering.numbering_scheme_transform[numbering.numbering_scheme][i];
			//TR << "Matched "<< landmark_to_match.chain() << " "<<landmark_to_match.resnum() << " " << landmark_to_match.insertion_code() << std::endl;
			//TR << "To " << new_landmark->chain() << " " << new_landmark->resnum() << " " << new_landmark->insertion_code()<< std::endl;
			return new_landmark;

		} else {
			continue;
		}
	}
	utility_exit_with_message("Could not get landmark. ");
}


void
AntibodyNumberingParser::read_cdr_definition_file(std::string const & file_path, AntibodyNumbering& numbering) {
	//Initializeour vector2d

	check_path(file_path);
	std::ifstream numbering_file(file_path.c_str());
	PyAssert((numbering_file.is_open()), "Unable to open CDR Definition file.");

	std::string line;

	for ( core::Size i = 1; i <= 6; ++i ) {
		vector1< PDBLandmarkOP > landmark_vec (CDRLandmarkEnum_total, nullptr);
		numbering.cdr_numbering.push_back(landmark_vec);
	}

	while ( getline(numbering_file, line) ) {

		//Skip any comments + empty lines
		if ( utility::startswith(line, "#") ) {
			continue;
		}
		if ( utility::startswith(line, "\n") ) {
			continue;
		}

		utility::trim(line, "\n"); //Remove trailing line break
		vector1< std::string > lineSP = utility::string_split_multi_delim(line); //Split on space or tab

		if ( lineSP[1] == "DEFINES" ) {
			read_cdr_definition_transform_line(lineSP, numbering);
		} else if ( enum_manager_->cdr_name_is_present(lineSP[1]) ) {
			read_cdr_definition_numbering_line(lineSP, numbering);
		} else {
			//TR << "Unrecognized CDR definition line.  Ignoring. " << line << std::endl;
			continue;
		}
	}
	TR << "Antibody CDR definition read successfully" << std::endl;
}

void
AntibodyNumberingParser::read_cdr_definition_transform_line(vector1<std::string> const & lineSP, AntibodyNumbering & numbering){
	for ( core::Size i = 2; i <= lineSP.size(); ++i ) {
		std::string definition_str = utility::string_split(lineSP[i], ':')[1];
		std::string scheme_used_str = utility::string_split(lineSP[i], ':')[2];
		if ( enum_manager_->cdr_definition_is_present(definition_str) && enum_manager_->numbering_scheme_is_present(scheme_used_str) ) {
			CDRDefinitionEnum cdr_definition = enum_manager_->cdr_definition_string_to_enum(definition_str);

			cdr_definitions_defined_.push_back(cdr_definition);
			cdr_definitions_defined_using_.push_back(enum_manager_->numbering_scheme_string_to_enum(scheme_used_str));
			//TR << definition_str << " "<< scheme_used_str << std::endl;
			for ( core::Size x = 1; x <= 6; ++x ) {
				vector1< PDBLandmarkOP > landmark_vec (CDRLandmarkEnum_total, nullptr);
				numbering.cdr_definition_transform[cdr_definition].push_back(landmark_vec);
			}
		} else {
			utility_exit_with_message("Unrecognized CDR Definition while parsing DEFINES line.  Please check: "+lineSP[i]);
		}
	}
	//TR << "AntibodyNumbering TRANSFORMS read successfully" << std::endl;
}

void
AntibodyNumberingParser::read_cdr_definition_numbering_line(vector1<std::string> const & lineSP, AntibodyNumbering & numbering) const{

	CDRNameEnum cdr = enum_manager_->cdr_name_string_to_enum(lineSP[1]);

	if ( ! enum_manager_->cdr_landmark_is_present(lineSP[2]) ) {
		utility_exit_with_message("Unrecognized antibody landmark: "+lineSP[2]);
	}
	CDRLandmarkEnum cdr_position = enum_manager_->cdr_landmark_string_to_enum(lineSP[2]);


	if ( cdr_definitions_defined_.size() + 2  !=  lineSP.size() ) {
		utility_exit_with_message("Transforms stated in DEFINES line but not defined in CDR line:"+lineSP[1]+" "+lineSP[2]);
	}

	for ( core::Size i = 1; i <= cdr_definitions_defined_.size(); ++i ) {

		core::Size lineSP_index = i+2;
		vector1< std::string > landmarkSP = utility::string_split(lineSP[lineSP_index], ':');


		char chain = landmarkSP[1].at(0);
		core::Size resnum; std::stringstream(landmarkSP[2]) >> resnum;
		char insertion_code;
		if ( landmarkSP[3] == "~" ) {
			insertion_code = ' ';
		} else {
			insertion_code = landmarkSP[3].at(0);
		}

		AntibodyNumberingSchemeEnum scheme = cdr_definitions_defined_using_[i];
		CDRDefinitionEnum definition = cdr_definitions_defined_[i];

		PDBLandmarkOP defined_pdb_landmark( new PDBLandmark(chain, resnum, insertion_code , scheme ) );
		//TR << "Access:" << defined_pdb_landmark->chain() << " " <<defined_pdb_landmark->resnum()<< " "<<defined_pdb_landmark->insertion_code()<< std::endl;
		PDBLandmarkOP new_landmark;


		//TR << enum_manager_->cdr_definition_enum_to_string(definition) << std::endl;

		if ( definition == numbering.cdr_definition && scheme == numbering.numbering_scheme ) {
			numbering.cdr_numbering[cdr][cdr_position] = defined_pdb_landmark;
			new_landmark = defined_pdb_landmark;
		} else if ( definition == numbering.cdr_definition ) {
			new_landmark = get_equivalent_landmark(numbering, defined_pdb_landmark->numbering_scheme(), *defined_pdb_landmark);
			numbering.cdr_numbering[cdr][cdr_position] = new_landmark;
		} else {
			new_landmark = get_equivalent_landmark(numbering, defined_pdb_landmark->numbering_scheme(), *defined_pdb_landmark);

		}

		numbering.cdr_definition_transform[definition][cdr][cdr_position] = new_landmark;

	}


}


/////PDBLandmark/////////
PDBLandmark::PDBLandmark(char chain, core::Size resnum, char insertion_code, AntibodyNumberingSchemeEnum scheme) {
	numbering_scheme_  =  scheme;
	resnum_ = resnum;
	chain_ = chain;
	insertion_code_ = insertion_code;
	//PDBLandmark(chain, resnum, insertion_code); This makes a temporary PDBLandmark that does nothing.  Blargh!

}

PDBLandmark::PDBLandmark(char chain, core::Size resnum, char insertion_code) {
	resnum_ = resnum;
	chain_ = chain;
	insertion_code_ = insertion_code;
}

PDBLandmark::~PDBLandmark()= default;


PDBLandmark &
PDBLandmark::operator =(const PDBLandmark& src){
	if ( this == &src ) {
		return *this;
	}
	numbering_scheme_ = src.numbering_scheme_;
	resnum_ = src.resnum_;
	chain_ = src.chain_;
	insertion_code_ = src.insertion_code_;

	return *this;
}

bool
PDBLandmark::operator ==(const PDBLandmark& compare) const {

	return (chain_ == compare.chain_ && resnum_ == compare.resnum_ && insertion_code_ == compare.insertion_code_);
}

bool
PDBLandmark::operator !=(const PDBLandmark& compare) const {
	return ! (chain_ == compare.chain_ && resnum_ == compare.resnum_ && insertion_code_ == compare.insertion_code_);
}

std::string
PDBLandmark::get_string() const {
	return utility::to_string(chain_)+":"+utility::to_string(resnum_)+":"+utility::to_string(insertion_code_);
}

} //antibody
} //protocols

