// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/metrics.cc
/// @brief Routines to measure antibodies
/// @author Jeffrey J. Gray

#ifndef INCLUDED_protocols_antibody_metrics_hh
#define INCLUDED_protocols_antibody_metrics_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <protocols/antibody/AntibodyInfo.hh>

namespace protocols {
namespace antibody {

using namespace core;
using namespace protocols::antibody;

/// @brief calculate the VH_VL packing angle from 2 sheet definitions on each antibody chain
utility::vector1< core::Real >
vl_vh_orientation_coords ( const core::pose::Pose & pose_in, const protocols::antibody::AntibodyInfo & ab_info );


///// kink measures /////

// @brief returns distance of the sc-sc Hbond across the strands at the beginning of the H3 kink (typically Asp-Arg)
core::Real
kink_RD_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);

// @brief returns distance for the bb-bb Hbond across the strands at the begining of the kink (typically Asp-Arg)
core::Real
kink_bb_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);


// @brief returns distance of the Trp sc-bb Hbond across the H3 kink residues (n-1 to n+2)
core::Real
kink_Trp_Hbond(const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo);

// @brief returns a pair of reals for the q distance and qbond dihedral angle from the four kink residues of the H3 C-terminal end
std::pair<core::Real,core::Real>
kink_dihedral( const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & abinfo, bool debug=false);


/// @brief calculate the SASA of the antibody paratope

std::pair<core::Real,core::Real>
paratope_sasa( const core::pose::Pose & pose, const protocols::antibody::AntibodyInfo & ab_info );

/// @brief calculate the net charge of the antibody
core::SSize
pose_charge( core::pose::Pose const & pose );

/// @brief calculate the net charge of the paratope
core::SSize
paratope_charge( core::pose::Pose const & pose, const protocols::antibody::AntibodyInfo & abinfo );

}
}


#endif

