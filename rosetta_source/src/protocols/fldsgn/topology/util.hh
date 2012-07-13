// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file ./src/protocols/src/fldsgn/topology/util.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_util_hh
#define INCLUDED_protocols_fldsgn_topology_util_hh

// unit header
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>
#include <protocols/fldsgn/topology/HelixPairing.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
#include <protocols/forge/build/Interval.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace topology {


/// @brief convert StrandParingSet of dssp to fldsgn::topology::StrandPairingSet
protocols::fldsgn::topology::StrandPairingSet
calc_strand_pairing_set(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::SS_Info2_COP const ssinfo,
	core::Size minimum_pair_length = 1
);


/// @brief calc delta sasa, when a molecule is splited to 2parts.
core::Real
calc_delta_sasa(
	core::pose::Pose const & pose,
	utility::vector1< protocols::forge::build::Interval > intervals,
	core::Real const pore_radius
);

/// @brief calculate chirality of three secondary structure elements
char
calc_chirality_sstriplet(
	SS_BaseCOP ss1,
	SS_BaseCOP ss2,
	SS_BaseCOP ss3
);
							 							 	
/// @brief check kink of helix, return number of loosen hydrogen
core::Size
check_kink_helix(
	 core::pose::Pose const & pose,
	 core::Size const begin,
	 core::Size const end 
);


/// @brief count number of internally made hydrogen bonds given a region
core::Size
count_internal_hbonds(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end
);

/// @brief count number of related hydrogen bonds given a region
core::Size
count_related_hbonds(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end
);
	
/// @brief get helix pairings
protocols::fldsgn::topology::HelixPairingSet	
calc_hpairset( 
   protocols::fldsgn::topology::SS_Info2_COP const ssinfo,
   core::Real distance = 15.0,
   core::Real angle = 45.0
);
	
} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
