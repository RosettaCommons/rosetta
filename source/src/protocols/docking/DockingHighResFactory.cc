// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/docking/DockingHighResFactory.cc
/// @brief Factory for creating DockingHighRes objects
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

// Unit headers
#include <protocols/docking/DockingHighResFactory.hh>

// Package headers
#include <protocols/docking/DockingHighRes.hh>
#include <protocols/moves/MoverFactory.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// C++ headers
#include <map>
#include <sstream>

namespace protocols {
namespace docking {

static THREAD_LOCAL basic::Tracer TR( "protocols.docking.DockingHighResFactory" );

/// @details Private constructor insures correctness of singleton.
DockingHighResFactory::DockingHighResFactory() {}

DockingHighResFactory::~DockingHighResFactory() = default;

/// @details The specified DHR_Type is mapped to the std::string that is known by the MoverFactory.
/// The MoverFactory generates a Mover instance, which is downcasted to a DockingHighResOP and returned.
DockingHighResOP
DockingHighResFactory::create_docking_high_res_mover(
	DHR_Type const & type
)
{
	using std::endl;
	using std::string;
	using std::stringstream;

	static const std::map<DHR_Type, string> type_map {
		{ DHR_Type::DockingHighResLegacy, "DockingHighResLegacy" },
		{ DHR_Type::DockingPrepackProtocol, "DockingPrepackProtocol" },
		{ DHR_Type::DockMCMProtocol, "DockMCMProtocol" },
		{ DHR_Type::DockMinMover, "DockMinMover" },
		{ DHR_Type::SnugDock, "SnugDock" },
		};

	// temporary until I decide how I want to handle this
	string type_name = type_map.at( type );

	TR.Trace << "generate DockingHighRes of type '" << type_name << "'." << endl;
	DockingHighResOP dhr( utility::pointer::dynamic_pointer_cast< DockingHighRes > ( moves::MoverFactory::get_instance()->newMover( type_name ) ) );
	if ( ! dhr ) {
		stringstream error_msg;
		error_msg
			<< "Attempting to create Mover "
			<< "'" << type_name << "' that is not a DockingHighRes." << endl
			<< "check spelling or "
			<< "register a new DockingHighRes in the MoverFactory" << endl;
		utility_exit_with_message( error_msg.str() );
	}
	return dhr;
}

} // namespace docking
} // namespace protocols
