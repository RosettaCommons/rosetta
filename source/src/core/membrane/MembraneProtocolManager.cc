// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MembraneProtocolManager.cc
///
/// @brief Protocol Config Manager
/// @detail Manage protocol-specific membrane changes and development
///
/// @author Rebecca Alford

#ifndef INCLUDED_core_membrane_MembraneProtocolManager_cc
#define INCLUDED_core_membrane_MembraneProtocolManager_cc

// Unit Headers
#include <core/membrane/MembraneProtocolManager.hh>

// Project Headers
#include <core/membrane/config/MembraneProtocolConfigLib.hh>

// Package Headers
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <cstdlib>
#include <string>

static basic::Tracer TR( "core.membrane.MembraneProtocolManager" );

namespace core {
namespace membrane {

/// @brief Empty Constructor
/// @param [none]
MembraneProtocolManager::MembraneProtocolManager() :
	utility::pointer::ReferenceCount(),
	is_managed_(false),
	is_ready_(false),
	protocol_(""),
	protocol_lib_(NULL)
{}

/// @brief Standard Constructor
/// @param protocol, mover
MembraneProtocolManager::MembraneProtocolManager(
	std::string protocol,
	protocols::moves::MoverOP mover
	) :
	utility::pointer::ReferenceCount(),
	is_managed_(false),
	is_ready_(false),
	protocol_(""),
	protocol_lib_(NULL)
{
	get_and_apply_changes( protocol, mover );
}

/// @brief Copy Constructor
/// @param Protocol manager object to copy to
MembraneProtocolManager::MembraneProtocolManager( MembraneProtocolManager const & src ) :
	utility::pointer::ReferenceCount()
{
	copy_data(*this, src);
}

/// @brief Destructor
/// @param [none]
MembraneProtocolManager::~MembraneProtocolManager()
{}

/// @brief Set is managed signal
/// @param is managed value
void MembraneProtocolManager::set_is_managed( bool is_managed ) { is_managed_ = is_managed; }

/// @brief Set ready signal
/// @param set ready value
void MembraneProtocolManager::set_is_ready( bool is_ready ) { is_ready_ = is_ready; }

/// @brief Set protocol type
/// @param protocol (str = legal values in option system)
void MembraneProtocolManager::set_protocol( std::string protocol ) { protocol_ = protocol; }

/// @brief Get is managed signal
/// @param [none]
bool MembraneProtocolManager::get_is_managed() { return is_managed_; }

/// @brief Get ready signal
/// @param [none]
bool MembraneProtocolManager::get_is_ready() { return is_ready_; }

/// @brief Get protocol name
/// @param [none]
std::string MembraneProtocolManager::get_protocol() { return protocol_; }

//////////////////////// Main Manager Helper Functions ////////////////
/// @brief Create protocol library function class
/// @param [none]
void MembraneProtocolManager::create_protocol_lib()
{
	TR << "Building protocol library" << std::endl;
	//protocol_lib_ = new MembraneProtocolConfigLib();
}

/// @brief Make Protocol Changes
/// @param protocol, mover
void MembraneProtocolManager::make_protocol_changes( std::string protocol, protocols::moves::MoverOP mover )
{
	TR << "Making protocol mover changes" << std::endl;
	/// Make standard protocol library function call (this probably wont work...)
	//protocol_lib_->protocol(mover);
	set_is_managed( true );
}

//////////////////////// Main Manager Helper Functions ////////////////
/// @brief Main protocol manager function
/// @param protocol, mover
void MembraneProtocolManager::get_and_apply_changes( std::string protocol, protocols::moves::MoverOP mover )
{
	TR << "Getting changes from protocol library and applying them to mover" << std::endl;

	/// Make Protocol library
	create_protocol_lib();

	// Make Protocol changes
	make_protocol_changes( protocol, mover );

	// Tell the membrane Hub you are done
	set_is_ready( true );

	TR << "Protocol is now managed by membrane protocol manager. Continuing..." << std::endl;

	// Done!
	return;
}

/// @brief Copy Constructor Helper Function
/// @param Object to copy from (lhs), object to copy to (rhs)
void MembraneProtocolManager::copy_data( MembraneProtocolManager lhs, MembraneProtocolManager rhs )
{
	// all of the copying
	lhs.is_managed_ = rhs.is_managed_;
	lhs.is_ready_ = rhs.is_ready_;

	lhs.protocol_ = rhs.protocol_;
	lhs.protocol_lib_ = rhs.protocol_lib_;
}

} // membrane
} // core

#endif INCLUDED_core_membrane_MembraneProtocolManager_cc