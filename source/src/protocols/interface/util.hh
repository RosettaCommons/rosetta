// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/interface/util.hh
/// @brief Collection of useful interface functions - both by chain and by groups. In one place.
///  The namespace is for interface classes and functions that are general-use and protocol level.
///  Many more specific classes and functions for interfaces can be found in their appropriate namespaces such as TaskOperations in protocols/toolbox/task_operations
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_interface_UTIL_HH
#define INCLUDED_protocols_interface_UTIL_HH

#include <core/pose/Pose.hh>

namespace protocols {
namespace interface {

//////////////// General Interface Tools ///////////////////////////

/// @brief Get a set of interface residues using the dock_chains interface string: Ex: LH_A.
/// @details Uses Steven's InterGroupNeighborsCalculator.  Does not require pose to have specific foldtree.
utility::vector1<bool>
select_interface_residues(core::pose::Pose const & pose, std::string interface, core::Size interface_distance);

//Add other general methods for selecting interface residues here using different metrics, etc.

//////////////// General Metrics ////////////////////////


}
}

#endif //#ifndef INCLUDED_protocols/interface_UTIL_HH

