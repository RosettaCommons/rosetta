// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file ./src/protocols/fldsgn/topology/Turn.fwd.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_Turn_fwd_hh
#define INCLUDED_protocols_fldsgn_topology_Turn_fwd_hh

#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

	class Turn;
	class TurnSet;

	typedef utility::pointer::owning_ptr< Turn > TurnOP;
	typedef utility::pointer::owning_ptr< Turn const > TurnCOP;
	typedef utility::vector1< TurnOP > Turns;
	typedef utility::pointer::owning_ptr< TurnSet > TurnSetOP;
	typedef utility::pointer::owning_ptr< TurnSet const > TurnSetCOP;


} // namespace topology
} // namespace fldsgn
} // namespace protocol

#endif
