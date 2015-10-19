// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

#ifndef INCLUDED_devel_sewing_util_HH
#define INCLUDED_devel_sewing_util_HH

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>

//Utility
#include <utility/vector1.hh>

//C++
#include <map>
#include <string>

namespace devel {
namespace sewing {

typedef std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > NativeRotamersMap;

//Create a file that saves native rotamers. Used later during design
void
dump_native_residue_file(
	NativeRotamersMap native_residue_map,
	std::string filename
);

} //sewing namespace
} //devel namespace

#endif
