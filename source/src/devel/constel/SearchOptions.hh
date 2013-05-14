// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Declaration of the different ways in which a pose is searched for
/// 	constellations.
/// @author Andrea Bazzoli

#ifndef INCLUDED_SearchOptions_hh
#define INCLUDED_SearchOptions_hh

#include <core/pose/Pose.fwd.hh>
#include <string>

using core::pose::Pose;

namespace devel {
namespace constel {

	/// @brief Search by target residue.
	void pair_constel_set(int const target_pdb_number, char const target_pdb_chain,
     Pose& pose_init);

	/// @brief Search by mutation pair.
	void pair_constel_set( std::string const& tgtmuts, Pose& pose_init );

} // constel
} // devel

#endif

