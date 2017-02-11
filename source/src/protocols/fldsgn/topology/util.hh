// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/src/fldsgn/topology/util.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_util_hh
#define INCLUDED_protocols_fldsgn_topology_util_hh

// unit header
#include <protocols/fldsgn/topology/StrandPairing.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
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


/// @brief check kink of helix, return number of loosen hydrogen
core::Size
check_kink_helix(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end );


utility::vector1< core::scoring::hbonds::HBond >
check_internal_hbonds(
	core::pose::Pose const & pose,
	core::Size const begin,
	core::Size const end );

core::Real
calc_strand_helix_angle(
	core::pose::Pose const & pose,
	protocols::fldsgn::topology::SS_Info2_COP const ssinfo,
	core::Size const strand_id1,
	core::Size const strand_id2,
	core::Size const helix_id,
	std::string const & geom_type
);

} // namespace topology
} // namespace fldsgn
} // namespace protocols

#endif
