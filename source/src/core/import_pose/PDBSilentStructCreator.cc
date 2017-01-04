// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/PDBSilentStructCreator.hh
/// @author James Thompson

// Unit Headers
#include <core/io/silent/SilentStructCreator.hh>

// Package Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/import_pose/PDBSilentStruct.hh>
#include <core/import_pose/PDBSilentStructCreator.hh>

#include <utility/vector1.hh>


namespace core {
namespace import_pose {

PDBSilentStructCreator::PDBSilentStructCreator() {}

PDBSilentStructCreator::~PDBSilentStructCreator() = default;

core::io::silent::SilentStructOP PDBSilentStructCreator::create_silent_struct( core::io::silent::SilentFileOptions const & opts ) const {
	return core::io::silent::SilentStructOP( new PDBSilentStruct( opts ) );
}

std::string PDBSilentStructCreator::keyname() const {
	return "pdb";
}

} //namespace import_pose
} //namespace core
