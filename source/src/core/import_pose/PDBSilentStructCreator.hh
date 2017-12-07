// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/BasicSilentStructCreator.hh
/// @brief  Base class for BasicSilentStructCreators for the BasicSilentStruct load-time factory registration scheme
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_PDBSilentStructCreator_hh
#define INCLUDED_core_import_pose_PDBSilentStructCreator_hh

// Unit Headers
#include <core/io/silent/SilentStructCreator.hh>

// c++ headers

#include <core/types.hh>


namespace core {
namespace import_pose {

/// @brief creator for the PDBSilentStruct class
class PDBSilentStructCreator : public core::io::silent::SilentStructCreator
{
public:
	PDBSilentStructCreator();
	~PDBSilentStructCreator() override;

	core::io::silent::SilentStructOP create_silent_struct( core::io::silent::SilentFileOptions const & opts ) const override;
	std::string keyname() const override;
};

} //namespace import_pose
} //namespace core

#endif
