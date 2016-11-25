// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/docking/DockingHighResFactory.hh
/// @brief Factory for creating DockingHighRes objects
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)


#ifndef INCLUDED_protocols_docking_DockingHighResFactory_HH
#define INCLUDED_protocols_docking_DockingHighResFactory_HH

// Unit headers
#include <protocols/docking/DockingHighResFactory.fwd.hh>
#include <utility/SingletonBase.hh>

// Package headers
#include <protocols/docking/DockingHighRes.fwd.hh>


namespace protocols {
namespace docking {

/// @brief These are the valid names of classes that can be instantiated by DockingHighResFactory
enum class DHR_Type {
	DockingHighResLegacy,
	DockingPrepackProtocol,
	DockMCMProtocol,
	DockMinMover,
	SnugDock
};

/// @brief DockingHighResFactory uses the MoverFactory to create instances of DockingHighRes
class DockingHighResFactory : public utility::SingletonBase< DockingHighResFactory >
{
public:
	friend class utility::SingletonBase< DockingHighResFactory >;

private:
	// Private constructor to make it singleton managed
	DockingHighResFactory();
	DockingHighResFactory(const DockingHighResFactory & src); // unimplemented

	DockingHighResFactory const &
	operator=( DockingHighResFactory const & ); // unimplemented

public:
	// Warning this is not called because of the singleton pattern
	virtual ~DockingHighResFactory();

	/// @brief Create the specified DockingHighRes
	DockingHighResOP
	create_docking_high_res_mover(
		DHR_Type const & type_name
	);
};

} // namespace docking
} // namespace protocols

#endif // INCLUDED_protocols_docking_DockingHighResFactory_HH
