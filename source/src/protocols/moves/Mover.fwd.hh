// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_Mover_fwd_hh
#define INCLUDED_protocols_moves_Mover_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

#include <map>

// Package headers

namespace protocols {
namespace moves {

class Mover;
typedef utility::pointer::shared_ptr< Mover > MoverOP;
typedef utility::pointer::shared_ptr< Mover const > MoverCOP;
typedef utility::pointer::weak_ptr< Mover > MoverAP;
typedef utility::pointer::weak_ptr< Mover const > MoverCAP;
typedef boost::shared_ptr< Mover> MoverSP;

typedef std::map< std::string const, MoverOP > Movers_map;

typedef std::map< std::string, std::string > SerializableState;
typedef boost::shared_ptr< SerializableState > SerializableStateSP;
typedef boost::shared_ptr< const SerializableState > SerializableStateCSP;
typedef boost::weak_ptr< SerializableState > SerializableStateWP;
typedef boost::weak_ptr< const SerializableState > SerializableStateCWP;

typedef std::map< std::string, std::string > MoverCache;
typedef boost::shared_ptr< MoverCache > MoverCacheSP;

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_Mover_fwd_HH
