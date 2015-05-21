// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/CDRSeqDesignOptions.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/string_util.hh>
#include <utility/py/PyAssert.hh>
#include <basic/Tracer.hh>
//#include <utility/io/izstream.hh>

#include <basic/database/open.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <iostream>
#include <fstream>
#include <cctype>
#include <utility/io/izstream.hh>

#include <boost/algorithm/string.hpp>

static thread_local basic::Tracer TR("protocols.antibody.design.CDRSeqDesignOptionsParser");


namespace protocols {
namespace antibody {
namespace design {

	using namespace boost;
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using utility::io::izstream;
	using std::string;
	using utility::vector1;

CDRSeqDesignOptions::CDRSeqDesignOptions():
	utility::pointer::ReferenceCount()
{

}

CDRSeqDesignOptions::CDRSeqDesignOptions(CDRNameEnum cdr):
	utility::pointer::ReferenceCount(),
	cdr_(cdr)
{

}

CDRSeqDesignOptions::CDRSeqDesignOptions(CDRSeqDesignOptions const & src):
	utility::pointer::ReferenceCount(src),
	cdr_(src.cdr_),
	design_(src.design_),
	design_strategy_(src.design_strategy_)
{

}

CDRSeqDesignOptions::~CDRSeqDesignOptions() {
}

void
CDRSeqDesignOptions::set_defaults() {
	design_ = true;
	design_strategy_ = seq_design_profiles;
}

void
CDRSeqDesignOptions::design_strategy(SeqDesignStrategyEnum strategy){
	design_strategy_ = strategy;
}

void
CDRSeqDesignOptions::design(bool design) {
	design_ = design;
}

void
CDRSeqDesignOptions::set_cdr(CDRNameEnum cdr) {
	cdr_ = cdr;
}

CDRSeqDesignOptionsOP
CDRSeqDesignOptions::clone() const {
	return CDRSeqDesignOptionsOP( new CDRSeqDesignOptions(*this) );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

CDRSeqDesignOptionsParser::CDRSeqDesignOptionsParser():
	utility::pointer::ReferenceCount(),
	default_and_user_(false)
{
	ab_manager_ = AntibodyEnumManagerOP( new AntibodyEnumManager() );
}

CDRSeqDesignOptionsParser::~CDRSeqDesignOptionsParser() {}

utility::vector1<CDRSeqDesignOptionsOP>
CDRSeqDesignOptionsParser::parse_default_and_user_options(std::string filename) {
	utility::vector1<CDRSeqDesignOptionsOP> antibody_options;
	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		antibody_options.push_back( parse_default_and_user_options( cdr, filename ) );
	}
	return antibody_options;
}

CDRSeqDesignOptionsOP
CDRSeqDesignOptionsParser::parse_default_and_user_options(CDRNameEnum cdr, std::string filename) {

	cdr_options_ = CDRSeqDesignOptionsOP( new CDRSeqDesignOptions(cdr) );
	std::string path = basic::options::option [basic::options::OptionKeys::antibody::design::base_instructions]();
	default_and_user_ = true;
	parse_options(cdr, path);
	parse_options(cdr, filename);
	default_and_user_ = false;
	return cdr_options_->clone();

}

utility::vector1<CDRSeqDesignOptionsOP>
CDRSeqDesignOptionsParser::parse_options(std::string filename) {
	utility::vector1<CDRSeqDesignOptionsOP> antibody_options;
	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		antibody_options.push_back( parse_options( cdr, filename ) );
	}
	return antibody_options;
}

CDRSeqDesignOptionsOP
CDRSeqDesignOptionsParser::parse_options(CDRNameEnum cdr, std::string path) {

	using namespace utility;
	using namespace std;

	if (default_and_user_){
		cdr_options_->set_cdr(cdr);
	}
	else{
		cdr_options_ = CDRSeqDesignOptionsOP( new CDRSeqDesignOptions(cdr) );
	}

	instructions_path_ = path;

	check_path();
	//This is straight from C++ tutorials.
	string line;
	izstream instruction_file(instructions_path_);
	if(instruction_file.bad()){
		utility_exit_with_message("Unable to open grafting instruction file.");
	}
	TR <<"Reading "<<path << " for "<< ab_manager_->cdr_name_enum_to_string(cdr) << std::endl;
	while (getline(instruction_file, line)){

		//Skip any comments + empty lines
		utility::trim(line, "\n"); //Remove trailing line break
		boost::algorithm::trim(line); //Remove any whitespace on either side of the string

		//Continue to next line on empty string, comment
		if (startswith(line, "#") || startswith(line, "\n") || line.empty()  ||  (line.find_first_not_of(' ') == std::string::npos) ){
			continue;
		}


		vector1< string > lineSP = string_split_multi_delim(line); //Split on space or tab
		check_line_len(lineSP, 2);
		//TR << utility::to_string(lineSP) <<std::endl;
		//Everything besides comments needs to have a CDR or ALL associated with it.
		std::string cdr_type = lineSP[1];
		boost::to_upper(cdr_type);

		std::string mode = lineSP[2];
		boost::to_upper(mode);

		if (cdr_type == "ALL"){
			parse_cdr_option(mode, lineSP);
		}
		else if (ab_manager_->cdr_name_is_present(cdr_type)){
			if (ab_manager_->cdr_name_string_to_enum(cdr_type) == cdr){
				parse_cdr_option(mode, lineSP);
			}
		}
		else {
			//If expansion to chains, frameworks, etc.  Do it here.
			//We may have separate parsers for framework or L2.5 or whatever.
			//If its not a CDR, just skip it for now so we can have
			TR << "Unrecognized CDR: "<<cdr_type <<" skipping...."<<std::endl;
			continue;
		}
	}
	instruction_file.close();
	TR << "Instructions read successfully" <<std::endl;
	return cdr_options_->clone();
}

void
CDRSeqDesignOptionsParser::check_path() {
	using namespace std;
	izstream check( instructions_path_, ifstream::in);
	if (check.good()){return;}
	else{
		ifstream check2((basic::database::full_name(instructions_path_, false)).c_str(), ifstream::in);
		if (check2.good()){
			instructions_path_ = basic::database::full_name(instructions_path_);
			return;
		}
		else{
			utility_exit_with_message("Instructions file path not good.  Please check path.");
		}
	}
}

void
CDRSeqDesignOptionsParser::parse_cdr_option(std::string const mode, vector1<string>& lineSP) {

	if (mode == "DESIGN" || mode == "SEQDESIGN" || mode == "SEQ_DESIGN" || mode == "SEQUENCEDESIGN" || mode == "SEQUENCE_DESIGN"){
		check_line_len(lineSP, 3);
		std::string adjective = lineSP[3];
		boost::to_upper(adjective);
		parse_cdr_design_option(adjective, lineSP);
	}
	else {
		parse_cdr_general_option(lineSP);
	}

}

void
CDRSeqDesignOptionsParser::check_line_len(const vector1<string> & lineSP, const core::Size len_check) const {
	if (lineSP.size() < len_check){
		utility_exit_with_message("Could not parse design instructions. Line not long enough: "+utility::to_string(len_check)+" "+utility::to_string(lineSP));
	}
}

void
CDRSeqDesignOptionsParser::parse_cdr_general_option(vector1<string> & lineSP) {

	check_line_len(lineSP, 2);
	std::string setting = lineSP[2];
	if (setting == "FIX"){
		cdr_options_->design(false);
	}
	else if (setting == "ALLOW") {
		cdr_options_->design(true);
	}
}

void
CDRSeqDesignOptionsParser::parse_cdr_design_option(std::string const adjective, vector1< string> & lineSP){

	using namespace utility;

	if (adjective=="FIX"){
		cdr_options_->design(false);
	}
	else if (adjective == "ALLOW") {
		cdr_options_->design(true);
	}
	else if(adjective=="PROFILES" || adjective == "PROFILE"){
		check_line_len(lineSP, 4);
		std::string option = lineSP[4];
		boost::to_upper(option);
		set_cdr_design_profile_option(option);
	}
	else{
		utility_exit_with_message("Could not parse ab design instruction.  Unknown option: "+adjective);
	}
}

void
CDRSeqDesignOptionsParser::set_cdr_design_profile_option(std::string const option) {

	if ( option == "CONSERVATIVE" ||  option=="CONSERVED"){
		cdr_options_->design_strategy(seq_design_conservative);
	}
	else if (option == "PROFILE"){
		cdr_options_->design_strategy(seq_design_profiles);
	}
	else if (option == "NONE" || option == "BASIC_DESIGN" || option == "BASIC"){
		cdr_options_->design_strategy(seq_design_basic);
	}
}

}
}
}
