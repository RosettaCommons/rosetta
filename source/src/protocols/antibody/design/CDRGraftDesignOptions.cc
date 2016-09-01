// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/CDRGraftDesignOptions.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
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

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.design.CDRGraftDesignOptions");


namespace protocols {
namespace antibody {
namespace design {

using namespace core;
//using namespace boost;
using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using utility::io::izstream;
using std::string;
using utility::vector1;

CDRGraftDesignOptions::CDRGraftDesignOptions():
	utility::pointer::ReferenceCount()
{
	set_defaults();
}

CDRGraftDesignOptions::CDRGraftDesignOptions(CDRNameEnum cdr):
	utility::pointer::ReferenceCount(),
	cdr_(cdr)
{
	set_defaults();
}

CDRGraftDesignOptions::CDRGraftDesignOptions(CDRGraftDesignOptions const & src):
	utility::pointer::ReferenceCount(src),
	cdr_(src.cdr_),
	design_(src.design_),
	mintype_(src.mintype_),
	min_neighbor_sc_(src.min_neighbor_sc_),
	min_sc_(src.min_sc_),
	min_rb_(src.min_rb_),
	neighbor_min_(src.neighbor_min_),
	cdr_weight_(src.cdr_weight_)
{

}



CDRGraftDesignOptions::~CDRGraftDesignOptions() {}

CDRGraftDesignOptionsOP
CDRGraftDesignOptions::clone() const {
	return CDRGraftDesignOptionsOP( new CDRGraftDesignOptions(*this) );
}

void
CDRGraftDesignOptions::set_defaults(){
	neighbor_min_.clear();
	design_ = true;
	min_sc_ = true;
	min_neighbor_sc_ = true;
	mintype_ = repack;
	min_rb_ = false;
	cdr_weight_ = 1.0;
}

void
CDRGraftDesignOptions::set_cdr(CDRNameEnum cdr) {
	cdr_ = cdr;
}

void
CDRGraftDesignOptions::design(bool design){
	design_ = design;
}

void
CDRGraftDesignOptions::weight(Real weight){
	cdr_weight_ = weight;
}

void
CDRGraftDesignOptions::min_neighbor_sc(bool min_neighbor_sc){
	min_neighbor_sc_ = min_neighbor_sc;
}

void
CDRGraftDesignOptions::min_sc(bool min_sc) {
	min_sc_ = min_sc;
}

void
CDRGraftDesignOptions::min_rb(bool min_rb){
	min_rb_ = min_rb;
}

void
CDRGraftDesignOptions::mintype(MinTypeEnum mintype){
	mintype_ = mintype;
}

void
CDRGraftDesignOptions::neighbor_min(utility::vector1<CDRNameEnum> neighbor_min){
	neighbor_min_ = neighbor_min;
}

void
CDRGraftDesignOptions::neighbor_min_add(CDRNameEnum cdr) {
	neighbor_min_.push_back(cdr);
}

void
CDRGraftDesignOptions::neighbor_min_clear(){
	neighbor_min_.clear();
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

CDRGraftDesignOptionsParser::CDRGraftDesignOptionsParser():
	utility::pointer::ReferenceCount(),
	default_and_user_(false)
{
	ab_manager_ = AntibodyEnumManagerOP( new AntibodyEnumManager() );
}

CDRGraftDesignOptionsParser::~CDRGraftDesignOptionsParser() {}

utility::vector1<CDRGraftDesignOptionsOP>
CDRGraftDesignOptionsParser::parse_default_and_user_options(std::string const & filename) {
	utility::vector1<CDRGraftDesignOptionsOP> antibody_options;
	for ( core::Size i = 1; i <= 6; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		antibody_options.push_back( parse_default_and_user_options( cdr, filename ) );
	}
	return antibody_options;
}

CDRGraftDesignOptionsOP
CDRGraftDesignOptionsParser::parse_default_and_user_options(CDRNameEnum cdr, std::string const & filename) {

	cdr_options_ = CDRGraftDesignOptionsOP( new CDRGraftDesignOptions(cdr) );
	std::string path = basic::options::option [basic::options::OptionKeys::antibody::design::base_cdr_instructions]();
	default_and_user_ = true;
	parse_options(cdr, path);
	parse_options(cdr, filename);
	default_and_user_ =false;
	return cdr_options_->clone();

}

utility::vector1<CDRGraftDesignOptionsOP>
CDRGraftDesignOptionsParser::parse_options(std::string const & filename) {
	utility::vector1<CDRGraftDesignOptionsOP> antibody_options;
	for ( core::Size i = 1; i <= 6; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		antibody_options.push_back( parse_options( cdr, filename ) );
	}
	return antibody_options;
}

CDRGraftDesignOptionsOP
CDRGraftDesignOptionsParser::parse_options(CDRNameEnum cdr, std::string const & path) {

	using namespace utility;
	using namespace std;

	if ( default_and_user_ ) {
		cdr_options_->set_cdr(cdr);
	} else {
		cdr_options_ = CDRGraftDesignOptionsOP( new CDRGraftDesignOptions(cdr) );
	}

	instructions_path_ = path;

	check_path();
	//This is straight from C++ tutorials.
	string line;
	izstream instruction_file(instructions_path_);
	if ( instruction_file.bad() ) {
		utility_exit_with_message("Unable to open grafting instruction file.");
	}
	//TR <<"Reading "<<path << " for "<< ab_manager_->cdr_name_enum_to_string(cdr) << std::endl;
	while ( getline(instruction_file, line) ) {

		//Skip any comments + empty lines
		utility::trim(line, "\n"); //Remove trailing line break
		boost::algorithm::trim(line); //Remove any whitespace on either side of the string

		//Continue to next line on empty string, comment
		if ( startswith(line, "#") || startswith(line, "\n") || line.empty()  ||  (line.find_first_not_of(' ') == std::string::npos) ) {
			continue;
		}

		boost::to_upper(line); //Capitalize everything.

		vector1< string > lineSP = string_split_multi_delim(line); //Split on space or tab
		check_line_len(lineSP, 2);
		//TR << utility::to_string(lineSP) << std::endl;

		//Everything besides comments needs to have a CDR or ALL associated with it.
		std::string cdr_type = lineSP[1];
		std::string mode = lineSP[2];


		if ( cdr_type == "ALL" && !(cdr == l4 || cdr == h4) ) {
			parse_cdr_option(mode, lineSP);
		} else if ( ( cdr_type == "DE" || cdr_type == "CDR4") && (cdr == l4 || cdr == h4) ) {
			parse_cdr_option(mode, lineSP);
		} else if ( ab_manager_->cdr_name_is_present(cdr_type) ) {
			if ( ab_manager_->cdr_name_string_to_enum(cdr_type) == cdr ) {
				parse_cdr_option(mode, lineSP);
			}
		} else {
			//If expansion to chains, frameworks, etc.  Do it here.
			//We may have separate parsers for framework or L2.5 or whatever.
			//If its not a CDR, just skip it for now so we can have
			//TR << "Unrecognized CDR: "<<cdr_type <<" skipping...."<<std::endl;
			continue;
		}
	}
	instruction_file.close();
	//TR << "Instructions read successfully" <<std::endl;

	if ( cdr == l4 || cdr == h4 ) {
		cdr_options_->design( false ); // Disable graft design of CDR4 as it is not yet implemented.
	}

	return cdr_options_->clone();
}

void
CDRGraftDesignOptionsParser::check_path() {
	using namespace std;
	izstream check( instructions_path_, ifstream::in);
	if ( check.good() ) { return;}
	else {
		ifstream check2((basic::database::full_name(instructions_path_, false)).c_str(), ifstream::in);
		if ( check2.good() ) {
			instructions_path_ = basic::database::full_name(instructions_path_);
			return;
		} else {
			utility_exit_with_message("Instructions file path not good.  Please check path.");
		}
	}
}

void
CDRGraftDesignOptionsParser::parse_cdr_option(std::string const & mode, vector1<string> const & lineSP) {

	///Since we are now including SeqDesign and GraftDesign as part of the overall protocol,
	/// it does not make sense to include MIN as part of GraftDesign, since it is overall part of the protocol.
	/// that said, I cannot make up my mind as to what the hell the the tag should be called.
	if ( mode == "GRAFT" || mode == "GRAFT_DESIGN" || mode == "GRAFTDESIGN" ||
			mode == "GRAFTING" || mode == "MINPROTOCOL" || mode == "MIN_PROTOCOL" || mode == "MIN_STEP" || mode == "MINSTEP" || mode == "MIN" ) {
		check_line_len(lineSP, 3);
		std::string adjective = lineSP[3];

		parse_cdr_graft_option(adjective, lineSP);
	} else if ( lineSP.size() == 2 ) {
		parse_cdr_general_option(lineSP);
	} else if ( lineSP.size() == 3 ) {
		std::string setting = lineSP[ 2 ];
		if ( setting == "WEIGHT" || setting == "WEIGHTS" ) {
			check_line_len(lineSP, 3);
			cdr_options_->weight(utility::string2Real(lineSP[3]));
		}
	}

}


void
CDRGraftDesignOptionsParser::check_line_len(vector1<string> const & lineSP, core::Size len_check) const {
	if ( lineSP.size() < len_check ) {
		utility_exit_with_message("Could not parse design instructions. Line not long enough: "+utility::to_string(len_check)+" "+utility::to_string(lineSP));
	}
}

void
CDRGraftDesignOptionsParser::parse_cdr_graft_option(std::string const & setting, vector1<string> const & lineSP) {

	//Here we match.  This is rather ugly, as I don't have much C++ expereince in this.  Python however...

	if ( lineSP.size() == 3 ) {
		std::string option = lineSP[3];
		set_cdr_graft_general_option(option);
		return;
	} else if ( setting == "MINTYPE" || setting == "MIN" ) {
		check_line_len(lineSP, 4);
		std::string type = lineSP[4];
		set_cdr_graft_mintype_options(type);
		return;
	} else if ( setting == "MIN_OTHER_CDRS" || setting == "MIN_NEIGHBOR_CDRS" || setting == "MIN_NEIGHBORS" || setting == "MIN_OTHER" || setting == "MINOTHER" ) {
		set_cdr_graft_neighbor_mintype_options(lineSP);
		return;
	} else if ( setting == "WEIGHT" || setting == "WEIGHTS" ) {
		check_line_len(lineSP, 4);
		cdr_options_->weight(utility::string2Real(lineSP[4]));
		return;
	} else {
		utility_exit_with_message("Cannot parse ab design instructions. " + setting);
	}
}

void
CDRGraftDesignOptionsParser::set_cdr_graft_general_option(std::string const & option) {

	if ( option == "FIX" ) {
		cdr_options_->design(false);
	} else if ( option == "ALLOW" ) {
		cdr_options_->design(true);
	} else { utility_exit_with_message("Unrecognized option: "+option); }
}

void
CDRGraftDesignOptionsParser::parse_cdr_general_option(vector1<string> const & lineSP) {

	check_line_len(lineSP, 2);
	std::string setting = lineSP[2];
	set_cdr_graft_general_option(setting);
}

void
CDRGraftDesignOptionsParser::set_cdr_graft_mintype_options(std::string const & mintype) {

	///If we ever have a design enum manager, add these to it to make this simpler.
	if ( mintype == "REPACK" || mintype == "PACK" ) {
		cdr_options_->mintype(protocols::antibody::design::repack);
	} else if ( mintype == "MIN" || mintype == "MINIMIZE" || mintype == "MINIMIZER" ) {
		cdr_options_->mintype(minimize);
	} else if ( mintype == "MIN_CART" || mintype == "MINIMIZE_CART" || mintype == "MINIMIZE_CARTESIAN" || mintype == "CARTMIN" ) {
		cdr_options_->mintype(minimize_cartesian);
	} else if ( mintype == "RELAX" ) {
		cdr_options_->mintype(relax);
	} else if ( mintype == "DUALSPACE" || mintype == "DUALSPACE_RELAX" || mintype == "DS_RELAX" ) {
		cdr_options_->mintype(dualspace);
	} else if ( mintype == "BACKRUB" ) {
		cdr_options_->mintype(backrub_protocol);
	} else if ( mintype == "INCLUDE_RB" || mintype == "RB" || mintype == "W_RB" ) {
		cdr_options_->min_rb(true);
	} else if ( mintype == "NONE" ) {
		cdr_options_->mintype(no_min);
	} else { utility_exit_with_message("Unrecognized mintype option: "+mintype); }

	//else if (mintype == "CENTROIDRELAX" || mintype == "CENTROID_RELAX" || mintype == "CENRELAX" || mintype == "CEN_RELAX"){
	// cdr_options_->mintype(centroid_relax);
	//}
	//else if (mintype == "CENTROIDRELAX_W_RB" || mintype == "CENTROID_RELAX_W_RB" || mintype=="CENRELAX_W_RB" || mintype == "CEN_RELAX_W_RB"){
	// graft_instructions_[cdr].mintype = centroid_relax;
	// graft_instructions_[cdr].min_rb = true;
	//}

}

void
CDRGraftDesignOptionsParser::set_cdr_graft_neighbor_mintype_options(utility::vector1<std::string> const & lineSP){

	///Clear neighbor min settings so overwrites do not get weird.
	cdr_options_->neighbor_min_clear();

	for ( core::Size i = 4; i <= lineSP.size(); ++i ) {
		std::string cdr_type = lineSP[i];
		if ( ! ab_manager_->cdr_name_is_present(cdr_type) ) {
			utility_exit_with_message("Unrecognized CDR type - "+cdr_type);
		}

		CDRNameEnum cdr = ab_manager_->cdr_name_string_to_enum(cdr_type);
		cdr_options_->neighbor_min_add(cdr);
	}
}


}
}
}











