// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterEnumManager.hh
/// @brief Functions for CDRClusterEnums
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_clusters_CDRCLUSTERENUMMANAGER_HH
#define INCLUDED_protocols_antibody_clusters_CDRCLUSTERENUMMANAGER_HH

#include <protocols/antibody/clusters/CDRClusterEnumManager.fwd.hh>

#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <map>
#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace antibody {
namespace clusters {

/// @brief Interface to this class is in AntibodyInfo.  Should be a singleton.
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

} //clusters
} //antibody
} //protocols

#endif //#ifndef INCLUDED_protocols/antibody_design_CDRCLUSTERENUMMANAGER_HH
