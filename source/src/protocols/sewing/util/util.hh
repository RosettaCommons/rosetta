// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocol_sewing_util_util_HH
#define INCLUDED_protocol_sewing_util_util_HH

//Package headers
#include <protocols/sewing/conformation/Assembly.hh>

//Protocol headers
#include <core/pose/Pose.hh>

//Utility
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/vector1.hh>

//Devel
#include <protocols/sewing/conformation/Model.fwd.hh>
#include <protocols/sewing/conformation/Assembly.fwd.hh>
#include <protocols/sewing/sampling/SewGraph.hh>

//C++
#include <map>
#include <string>

namespace protocols {
namespace sewing  {

std::map<core::id::AtomID, core::id::AtomID>
largest_continuous_atom_map(
	std::map<core::id::AtomID, core::id::AtomID> const & atom_map
);

} //sewing namespace
} //protocols namespace

#endif
