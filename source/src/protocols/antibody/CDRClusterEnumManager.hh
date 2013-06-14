// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/CDRClusterEnumManager.hh
/// @brief Functions for CDRClusterEnums
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_CDRCLUSTERENUMMANAGER_HH
#define INCLUDED_protocols_antibody_CDRCLUSTERENUMMANAGER_HH

#include <protocols/antibody/CDRClusterEnumManager.fwd.hh>

#include <protocols/antibody/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <map>
#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace antibody {


///@brief Interface to this class is in AntibodyInfo.
class CDRClusterEnumManager : public utility::pointer::ReferenceCount {
public:

	CDRClusterEnumManager();

	virtual ~CDRClusterEnumManager();

	CDRClusterEnum
	cdr_cluster_string_to_enum(std::string const & cluster) const;

	std::string
	cdr_cluster_enum_to_string(CDRClusterEnum const cluster) const;

	bool
	cdr_cluster_is_present(std::string const & cluster) const;

private:


	void setup();

	//bool initialized_;
	utility::vector1< std::string > enum_to_string_;
	std::map< std::string, CDRClusterEnum > string_to_enum_;

};


}
}

#endif	//#ifndef INCLUDED_protocols/antibody_design_CDRCLUSTERENUMMANAGER_HH
