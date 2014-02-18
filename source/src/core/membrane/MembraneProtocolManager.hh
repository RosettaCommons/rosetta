// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MembraneProtocolManager.hh
///
/// @brief Manage protocol specific changes
/// @detail Communicates with the membrane protocol config library to make changes to an apply in mover
///
/// @author Rebecca Alford

#ifndef INCLUDED_core_membrane_MembraneProtocolManager_hh
#define INCLUDED_core_membrane_MembraneProtocolManager_hh

// Unit Headers
#include <core/membrane/MembraneProtocolManager.fwd.hh>

// Project Headers
#include <core/membrane/config/MembraneProtocolConfigLib.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <algorithm>

namespace core {
namespace membrane {

/// @briec Membrane Protocol Manager Class
class MembraneProtocolManager : public utility::pointer::ReferenceCount {

public:

	/// @brief Empty Constructor
	/// @param [none]
	MembraneProtocolManager();

	/// @brief Standard Constructor
	/// @param protocol, mover
	MembraneProtocolManager(
		std::string protocol,
		protocols::moves::MoverOP mover
		);

	/// @brief Copy Constructor
	/// @param Protocol manager object to copy to
	MembraneProtocolManager( MembraneProtocolManager const & src );

	/// @brief Destructor
	/// @param [none]
	~MembraneProtocolManager();

	/// @brief Set is managed signal
	/// @param is managed value
	void set_is_managed( bool is_managed );

	/// @brief Set ready signal
	/// @param set ready value
	void set_is_ready( bool is_ready );

	/// @brief Set protocol type
	/// @param protocol (str = legal values in option system)
	void set_protocol( std::string protocol );

	/// @brief Get is managed signal
	/// @param [none]
	bool get_is_managed();

	/// @brief Get ready signal
	/// @param [none]
	bool get_is_ready();

	/// @brief Get protocol name
	/// @param [none]
	std::string get_protocol();

private:

	/// @brief Setup State Vars
	bool is_managed_;
	bool is_ready_;

	/// @brief Protocol Options
	std::string protocol_;

	/// @brief Protocol Change Library
	core::membrane::config::MembraneProtocolConfigLibOP protocol_lib_;

	/// @brief Create protocol library function class
	/// @param [none]
	void create_protocol_lib();

	/// @brief Make Protocol Changes
	/// @param protocol, mover
	void make_protocol_changes( std::string protocol, protocols::moves::MoverOP mover );

	/// @brief Main Protocol Manager funciton
	/// @param protocol, mover
	void get_and_apply_changes( std::string protocol, protocols::moves::MoverOP mover );

	/// @brief Copy Constructor Helper Function
	/// @param Object to copy from (lhs), object to copy to (rhs)
	void copy_data( MembraneProtocolManager lhs, MembraneProtocolManager rhs );

}; // class MembraneProtocolManager

} // membrane
} // core

#endif // INCLUDED_core_membrane_MembraneProtocolManager_hh