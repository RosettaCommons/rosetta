// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/SearchOptions.hh
/// @brief Declaration of the different ways in which a pose is searched for
///  constellations.
/// @author Andrea Bazzoli

#ifndef INCLUDED_SearchOptions_hh
#define INCLUDED_SearchOptions_hh

#include <core/pose/Pose.fwd.hh>
#include <string>

namespace protocols {
namespace constel {

class NeighTeller;

/// @brief Searches pair-constellations by target residue.
void pair_constel_set(int const target_pdb_number, std::string const & target_pdb_chain,
	core::pose::Pose& pose_init);

/// @brief Searches pair-constellations by mutation pair.
void pair_constel_set( std::string const& tgtmuts, core::pose::Pose& pose_init );

/// @brief Searches triple-constellations by target residue.
void triple_constel_set(int const target_pdb_number, std::string const & target_pdb_chain,
	core::pose::Pose& pose_init);

/// @brief Searches a single, target constellation.
void target_constel(std::string &tgtcnl_fil, core::pose::Pose& ps);

} // constel
} // protocols

#endif

