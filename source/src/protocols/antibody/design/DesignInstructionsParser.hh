// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/DesignInstructionsParser.hh
/// @brief Parsers Antibody Design Instruction files
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_DESIGNINSTRUCTIONSPARSER_HH
#define INCLUDED_protocols_antibody_design_DESIGNINSTRUCTIONSPARSER_HH

#include <protocols/antibody/design/AntibodyGraftDesignMover.hh>
#include <protocols/antibody/design/AntibodySeqDesignMover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace antibody {
namespace design {

	using namespace utility;
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using std::string;
	
///@brief Lightweight class for parsing design instructions.
class DesignInstructionsParser{

	
typedef std::map< CDRNameEnum, CDRGraftInstructions > GraftInstructions;
typedef std::map< CDRNameEnum, CDRDesignInstructions > DesignInstructions;

public:
	DesignInstructionsParser(AntibodyInfoCOP ab_info, string const path);
	
	~DesignInstructionsParser();

	void
	read_cdr_graft_instructions(GraftInstructions & instructions);
	
	void
	read_cdr_design_instructions(DesignInstructions & instructions);
	
private:
	
	///Tries to find the path in either database, relative, or absolute.
	void
	check_path() ;
	
	///Graft Instructions Methods
	//void
	//parse_cdr_graft_weights(GraftInstructions & instructions, vector1< string > & lineSP);
	
	void
	parse_cdr_graft_general_options(GraftInstructions & instructions, vector1< string > & lineSP) const;
	
	void
	parse_cdr_graft_type_options(GraftInstructions & instructions, vector1 < string > & lineSP) const;
	
	void
	parse_cdr_graft_mintype(GraftInstructions & instructions, vector1 <string> & lineSP) const;
	
	void
	parse_cdr_design_option(DesignInstructions & instructions, vector1 <string> & lineSP) const;
	
	std::string instructions_path_;
	AntibodyInfoCOP ab_info_;
	AntibodyEnumManagerOP ab_manager_;
	CDRClusterEnumManagerOP cluster_manager_;
}; 
}
}
}

#endif //INCLUDED_protocols/antibody/design_DESIGNINSTRUCTIONSPARSER_HH


