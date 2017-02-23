// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyDesignEnumManager.hh
/// @brief Functions for AntibodyEnumerators
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignEnumManager_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignEnumManager_hh

#include <protocols/antibody/design/AntibodyDesignEnumManager.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <utility/vector1.hh>
#include <map>
#include <string>

namespace protocols {
namespace antibody {
namespace design {


class AntibodyDesignEnumManager : public utility::pointer::ReferenceCount {


public:

	AntibodyDesignEnumManager();

	virtual ~AntibodyDesignEnumManager();
	
	////////////////// SeqDesignStrategy ///////////////////////////////////////////////

	SeqDesignStrategyEnum
	seq_design_strategy_string_to_enum(std::string const &seq_design_strategy) const;

	std::string
	seq_design_strategy_enum_to_string(SeqDesignStrategyEnum const seq_design_strategy) const;

	bool
	seq_design_strategy_is_present(std::string const &seq_design_strategy) const;



	////////////////// Design Protocol ///////////////////////////////////////////////

	AntibodyDesignProtocolEnum
	design_protocol_string_to_enum(std::string const &design_protocol) const;

	std::string
	design_protocol_enum_to_string(AntibodyDesignProtocolEnum const design_protocol) const;

	bool
	design_protocol_is_present(std::string const &design_protocol) const;

	utility::vector1<std::string>
	get_recognized_design_protocols() const;


	////////////////// MinType ///////////////////////////////////////////////

	MinTypeEnum
	min_type_string_to_enum(std::string const &min_type) const;

	std::string
	min_type_enum_to_string(MinTypeEnum const min_type) const;

	bool
	min_type_is_present(std::string const &min_type) const;


private:


	void setup();


	utility::vector1<std::string> seq_design_strategy_to_string_;
	std::map<std::string, SeqDesignStrategyEnum> seq_design_strategy_to_enum_;

	utility::vector1<std::string> design_protocol_to_string_;
	std::map<std::string, AntibodyDesignProtocolEnum> design_protocol_to_enum_;

	utility::vector1<std::string> min_type_to_string_;
	std::map<std::string, MinTypeEnum> min_type_to_enum_;


};

}
}
}

#endif //#ifndef INCLUDED_protocols/antibody/ANTIBODYENUMMANAGER_HH

