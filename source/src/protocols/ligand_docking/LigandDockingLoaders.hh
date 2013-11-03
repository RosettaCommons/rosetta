// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandDockingLoader.hh
/// @brief  Declartion of the InterfaceBuilderLoader and MoveMapBuilderLoader classes
///         for adding named InterfaceBuilders and MoveMapBuilders into the parser's basic::datacache::DataMap
/// @author Gordon Lemmon
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) -- moved here from DockDesignParser.cc

#ifndef INCLUDED_protocols_ligand_docking_LigandDockingLoaders_hh
#define INCLUDED_protocols_ligand_docking_LigandDockingLoaders_hh

// Package Headers
#include <protocols/jd2/parser/DataLoader.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace ligand_docking {

/// @brief A class for loading InterfaceBuilders into the XML parser's basic::datacache::DataMap.
class InterfaceBuilderLoader : public jd2::parser::DataLoader
{
public:
	InterfaceBuilderLoader();
	virtual ~InterfaceBuilderLoader();

	/// @brief The InterfaceBuilderLoader will create named InterfaceBuilders and load them into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

};

/// @brief A class for loading MoveMapBuilders into the XML parser's basic::datacache::DataMap.
/// NOTE that in the input .xml file, the InterfaceBuilder must be specified before
/// the MoveMapBuilder
class MoveMapBuilderLoader : public jd2::parser::DataLoader
{
public:
	MoveMapBuilderLoader();
	virtual ~MoveMapBuilderLoader();

	/// @brief The InterfaceBuilderLoader will create named InterfaceBuilders and load them into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

};

class LigandAreaLoader : public jd2::parser::DataLoader
{
public:
	LigandAreaLoader();
	virtual ~LigandAreaLoader();

	/// @brief The InterfaceBuilderLoader will create named InterfaceBuilders and load them into the basic::datacache::DataMap
	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) const;

};

} //namespace ligand_docking
} //namespace protocols

#endif
