// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyDesignEnumManager.cc
/// @brief Functions for AntibodyEnumerators
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyDesignEnumManager.hh>

#include <map>
#include <string>

#include <utility/vector1.hh>
#include <utility/py/PyAssert.hh>
#include <utility/string_util.hh>
#include <boost/algorithm/string.hpp>

namespace protocols {
namespace antibody {
namespace design {


AntibodyDesignEnumManager::AntibodyDesignEnumManager() {
	setup();
}

AntibodyDesignEnumManager::~AntibodyDesignEnumManager() { }

void
AntibodyDesignEnumManager::setup() {

	///Initialize the length of the vectors
	seq_design_strategy_to_string_.resize(SeqDesignStrategyEnum_total);
	min_type_to_string_.resize(MinTypeEnum_total);
	design_protocol_to_string_.resize(DesignProtocolEnum_total);

	///Manually construct the variables

	////////////////// Seq Design Strategy ///////////////////////////////////////////////
	seq_design_strategy_to_string_[seq_design_conservative] = "seq_design_conservative";
	seq_design_strategy_to_string_[seq_design_profiles] = "seq_design_profiles";
	seq_design_strategy_to_string_[seq_design_basic] = "seq_design_basic";
	seq_design_strategy_to_string_[seq_design_none] = "seq_design_none";


	seq_design_strategy_to_enum_["CONSERVATIVE"] = seq_design_conservative;
	seq_design_strategy_to_enum_["CONSERVED"] = seq_design_conservative;
	seq_design_strategy_to_enum_["SEQ_DESIGN_CONSERVATIVE\""] = seq_design_conservative;

	seq_design_strategy_to_enum_["PROFILE"] = seq_design_profiles;
	seq_design_strategy_to_enum_["PROFILES"] = seq_design_profiles;
	seq_design_strategy_to_enum_["SEQ_DESIGN_PROFILES"] = seq_design_profiles;

	seq_design_strategy_to_enum_["PROFILE_SETS"] = seq_design_profile_sets;
	seq_design_strategy_to_enum_["PROFILESETS"] = seq_design_profile_sets;
	seq_design_strategy_to_enum_["SEQ_DESIGN_PROFILE_SETS"] = seq_design_profile_sets;

	seq_design_strategy_to_enum_["PROFILES_COMBINED"] = seq_design_profile_sets_combined;
	seq_design_strategy_to_enum_["PROFILES_AND_SETS"] = seq_design_profile_sets_combined;
	seq_design_strategy_to_enum_["PROFILES_AND_PROFILE_SETS"] = seq_design_profile_sets_combined;

	seq_design_strategy_to_enum_["BASIC_DESIGN"] = seq_design_basic;
	seq_design_strategy_to_enum_["BASIC"] = seq_design_basic;
	seq_design_strategy_to_enum_["NORMAL"] = seq_design_basic;
	seq_design_strategy_to_enum_["SEQ_DESIGN_BASIC"] = seq_design_basic;

	seq_design_strategy_to_enum_["SEQ_DESIGN_NONE"] = seq_design_none;
	seq_design_strategy_to_enum_["NO_DESIGN"] = seq_design_none;
	seq_design_strategy_to_enum_["OFF"] = seq_design_none;


	////////////////// Design Protocol ///////////////////////////////////////////////
	design_protocol_to_string_[generalized_monte_carlo] = "generalized_monte_carlo";
	design_protocol_to_string_[deterministic_graft] = "deterministic_graft";
	design_protocol_to_string_[even_cluster_monte_carlo] = "even_cluster_monte_carlo";
	design_protocol_to_string_[even_length_cluster_monte_carlo] = "even_length_cluster_monte_carlo";

	design_protocol_to_enum_["GENERALIZED_MONTE_CARLO"] = generalized_monte_carlo;
	design_protocol_to_enum_["GEN_MC"] = generalized_monte_carlo;

	design_protocol_to_enum_["DETERMINISTIC_GRAFT"] = deterministic_graft;

	design_protocol_to_enum_["EVEN_CLUSTER_MC"] = even_cluster_monte_carlo;

	design_protocol_to_enum_["EVEN_LENGTH_CLUSTER_MC"] = even_length_cluster_monte_carlo;


	////////////////// MinType ///////////////////////////////////////////////
	min_type_to_string_[relax] = "relax";
	min_type_to_string_[centroid_relax] = "centroid_relax";
	min_type_to_string_[minimize] = "minimize";
	min_type_to_string_[minimize_cartesian] = "minimize_cartesian";
	min_type_to_string_[dualspace] = "dualspace";
	min_type_to_string_[repack] = "repack";
	min_type_to_string_[backrub_protocol] = "backrub_protocol";
	min_type_to_string_[no_min] = "no_min";

	min_type_to_enum_["RELAX"] = relax;
	min_type_to_enum_["CENTROID_RELAX"] = centroid_relax;
	min_type_to_enum_["MINIMIZE"] = minimize;
	min_type_to_enum_["MINIMIZE_CARTESIAN"] = minimize_cartesian;
	min_type_to_enum_["DUALSPACE"] = dualspace;
	min_type_to_enum_["REPACK"] = repack;
	min_type_to_enum_["BACKRUB_PROTOCOL"] = backrub_protocol;
	min_type_to_enum_["NO_MIN"] = no_min;


}

SeqDesignStrategyEnum
AntibodyDesignEnumManager::seq_design_strategy_string_to_enum(std::string const &seq_design_strategy) const {

	std::string s = boost::to_upper_copy<std::string>(seq_design_strategy);

	if ( ! seq_design_strategy_is_present( s ) ) {
		utility_exit_with_message("Strategy not present!" +s);
	}

	std::map<std::string, SeqDesignStrategyEnum>::const_iterator iter(
		seq_design_strategy_to_enum_.find(s));
	return iter->second;
}

std::string
AntibodyDesignEnumManager::seq_design_strategy_enum_to_string(SeqDesignStrategyEnum const seq_design_strategy) const {
	return seq_design_strategy_to_string_[seq_design_strategy];
}

bool
AntibodyDesignEnumManager::seq_design_strategy_is_present(std::string const &seq_design_strategy) const {

	std::map<std::string, SeqDesignStrategyEnum>::const_iterator iter(
		seq_design_strategy_to_enum_.find(seq_design_strategy));
	return iter != seq_design_strategy_to_enum_.end();
}


AntibodyDesignProtocolEnum
AntibodyDesignEnumManager::design_protocol_string_to_enum(std::string const &design_protocol) const {

	std::string s = boost::to_upper_copy<std::string>(design_protocol);
	if ( ! design_protocol_is_present( s ) ) {
		utility_exit_with_message("Design Protocol not present! "+ s );
	}

	std::map<std::string, AntibodyDesignProtocolEnum>::const_iterator iter(design_protocol_to_enum_.find(s));
	return iter->second;
}

std::string
AntibodyDesignEnumManager::design_protocol_enum_to_string(AntibodyDesignProtocolEnum const design_protocol) const {
	return design_protocol_to_string_[design_protocol];
}

bool
AntibodyDesignEnumManager::design_protocol_is_present(std::string const & design_protocol) const {

	std::map<std::string, AntibodyDesignProtocolEnum>::const_iterator iter(
		design_protocol_to_enum_.find(design_protocol));
	return iter != design_protocol_to_enum_.end();
}


MinTypeEnum
AntibodyDesignEnumManager::min_type_string_to_enum(std::string const & min_type) const {

	std::string s = boost::to_upper_copy<std::string>(min_type);

	if ( ! min_type_is_present( s ) ) {
		utility_exit_with_message("MinType not present! "+s);
	}
	std::map<std::string, MinTypeEnum>::const_iterator iter(min_type_to_enum_.find(s));
	return iter->second;
}

std::string
AntibodyDesignEnumManager::min_type_enum_to_string(MinTypeEnum const min_type) const {
	return min_type_to_string_[min_type];
}

bool
AntibodyDesignEnumManager::min_type_is_present(std::string const & min_type) const {

	std::map<std::string, MinTypeEnum>::const_iterator iter(min_type_to_enum_.find(min_type));
	return iter != min_type_to_enum_.end();
}

}
}
}

