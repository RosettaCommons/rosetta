// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/CDRSetOptionsParser.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_database_CDRSetOptionsParser_hh
#define INCLUDED_protocols_antibody_database_CDRSetOptionsParser_hh

#include <protocols/antibody/database/CDRSetOptionsParser.fwd.hh>

#include <protocols/antibody/database/CDRSetOptions.fwd.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.fwd.hh>
#include <protocols/antibody/AntibodyEnumManager.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
			

///@brief Parses CDRSetOptions for a single CDR at a time from a simple instruction file
class CDRSetOptionsParser : public utility::pointer::ReferenceCount {
public:

	CDRSetOptionsParser();

	virtual ~CDRSetOptionsParser();
	
	CDRSetOptionsOP
	parse_options(CDRNameEnum cdr, std::string filename);
	
	///@brief Parse default instructions (mainly used for AbDesign) then parse user file
	CDRSetOptionsOP
	parse_default_and_user_options(CDRNameEnum cdr, std::string filename);
	
	//CDRSetOptionsOP
	//parse_default_and_user_options(CDRNameEnum cdr);
	
	
	/////////////Parse Options for all CDRs////////////////
	
	utility::vector1<CDRSetOptionsOP>
	parse_options(std::string filename);
	
	utility::vector1<CDRSetOptionsOP>
	parse_default_and_user_options(std::string filename);

	//utility::vector1<CDRSetOptionsOP>
	//parse_default_and_user_options();
	
	
private:
	
	///@brief Tries to find the path in either database, relative, or absolute, sets it.
	void
	check_path();
	
	void
	check_line_len(const utility::vector1<std::string> & lineSP, const core::Size len_check) const;
	
	void
	parse_cdr_option(std::string const mode, utility::vector1< std::string> & lineSP);
	
	void
	parse_cdr_set_option(std::string const adjective, utility::vector1<std::string> & lineSP) ;
	
	void
	set_cdr_set_exclude_options(std::string const type, utility::vector1<std::string> & lineSP);
	
	void
	clear_cdr_set_exclude_options(std::string const type);
	
	void
	set_cdr_set_include_options(std::string const type, utility::vector1<std::string> & lineSP);
	
	void
	clear_cdr_set_include_options(std::string const type);
	
	void
	set_cdr_set_mintype_options(std::string const mintype);
	
	void
	set_cdr_set_general_option(std::string const option);
	
private:
	
	std::string instructions_path_;
	CDRSetOptionsOP cdr_options_;
	AntibodyEnumManagerCOP ab_manager_;
	clusters::CDRClusterEnumManagerCOP cluster_manager_;
	bool default_and_user_;
};

}
}




#endif	





