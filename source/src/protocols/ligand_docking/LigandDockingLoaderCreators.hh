// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandDockingLoaderCreators.hh
/// @brief  Creator classes for the ligand docking DataLoader classes: InterfaceBuilderLoader and MoveMapBuilderLoader
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_ligand_docking_LigandDockingLoaderCreators_hh
#define INCLUDED_protocols_ligand_docking_LigandDockingLoaderCreators_hh

// Package headers
#include <protocols/jd2/parser/DataLoaderCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace ligand_docking {

class InterfaceBuilderLoaderCreator : public jd2::parser::DataLoaderCreator
{
public:
	jd2::parser::DataLoaderOP create_loader() const override;
	std::string keyname() const override;
};

class MoveMapBuilderLoaderCreator : public jd2::parser::DataLoaderCreator
{
public:
	jd2::parser::DataLoaderOP create_loader() const override;
	std::string keyname() const override;
};

class LigandAreaLoaderCreator : public jd2::parser::DataLoaderCreator
{
public:
	jd2::parser::DataLoaderOP create_loader() const override;
	std::string keyname() const override;
};

} //namespace ligand_docking
} //namespace protocols

#endif
