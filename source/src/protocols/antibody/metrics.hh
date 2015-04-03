// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/metrics.cc
/// @brief Routines to measure antibodies
/// @author Jeffrey J. Gray, Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_antibody_metrics_hh
#define INCLUDED_protocols_antibody_metrics_hh

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

namespace protocols {
namespace antibody {

using namespace core;

/// @brief calculate the VH_VL packing angle from 2 sheet definitions on each antibody chain
/// @details vector is made up of:
///  vl_vh_distance
///  opening_angle
///  opposite_angle
///  packing_angle
utility::vector1< core::Real >
vl_vh_orientation_coords ( const core::pose::Pose & pose_in, const AntibodyInfo & ab_info );


///// kink measures /////

// @brief returns distance of the sc-sc Hbond across the strands at the beginning of the H3 kink (typically Asp-Arg)
core::Real
kink_RD_Hbond(const core::pose::Pose & pose, const AntibodyInfo & abinfo);

// @brief returns distance for the bb-bb Hbond across the strands at the begining of the kink (typically Asp-Arg)
core::Real
kink_bb_Hbond(const core::pose::Pose & pose, const AntibodyInfo & abinfo);


// @brief returns distance of the Trp sc-bb Hbond across the H3 kink residues (n-1 to n+2)
core::Real
kink_Trp_Hbond(const core::pose::Pose & pose, const AntibodyInfo & abinfo);

// @brief returns a pair of reals for the q distance and qbond dihedral angle from the four kink residues of the H3 C-terminal end
std::pair<core::Real,core::Real>
kink_dihedral( const core::pose::Pose & pose, const AntibodyInfo & abinfo, bool debug=false);


///// paratope and pose measures /////

/// @brief Convenience struct for passing around paratope data, including individual data for cdrs.
/// CDRs not present default to 0.  Templates are used mainly to have Size/ Real or string values. 
template <typename T>
struct ParatopeMetric {

    ParatopeMetric():
        paratope(0),
        paratope_heavy(0),
        paratope_light(0)
    {}

	vector1< T > cdr;
	T paratope;
	T paratope_heavy;
	T paratope_light;
};

/// @brief calculate the SASA of the antibody paratope.  Returns all and hydrophobic components.
std::pair<ParatopeMetric <core::Real>, ParatopeMetric<core::Real> >
paratope_sasa( const core::pose::Pose & pose, const AntibodyInfo & ab_info);

/// @brief calculate the net charge of the paratope
ParatopeMetric<core::SSize>
paratope_charge( core::pose::Pose const & pose, const protocols::antibody::AntibodyInfo & abinfo );

 /// @brief calculate the net charge of the antibody
core::SSize
pose_charge( core::pose::Pose const & pose );

/// @brief calculate dSASA of the paratope and antigen using result of InterfaceAnalyzerMover
//std::pair<core::Real,core::Real>
//paratope_dSASA( 
//	const core::pose::Pose & pose,
//	const AntibodyInfo & ab_info);

/// @brief calculate dSASA of the paratope using data held in PerResidueInterfaceData returned from InterfaceAnalyzerMover
//std::pair<core::Real,core::Real>
//paratope_dSASA( 
//	const core::pose::Pose & pose,
//	const AntibodyInfo & ab_info
//	const protocols::analysis::PerResidueInterfaceData & interface_data);

/// @brief calculate dSASA of a cdr and antigen using InterfaceAnalyzerMover. Returns all and hydrophobic components.
//std::pair<core::Real, core::Real>
//cdr_dSASA(
//	const core::pose::Pose & pose,
//	const AntibodyInfo & ab_info,
//	const CDRNameEnum & cdr);

/// @brief calculate dSASA of a cdr using data held in PerResidueInterfaceData returned from InterfaceAnalyzerMover.
//std::pair<core::Real, core::Real>
//cdr_dSASA(
//	const core::pose::Pose & pose,
//	const AntibodyInfo & ab_info,
//	const CDRNameEnum & cdr,
//	const protocols::analysis::PerResidueInterfaceData);

///// CDR measures /////

/// @brief Calculates energy of cdr by decomposing scorefxn into pair and then summing.  This means it includes hbonding components, etc.
core::Real
cdr_energy(core::pose::Pose const & pose, AntibodyInfoCOP ab_info, core::scoring::ScoreFunctionCOP scorefxn, CDRNameEnum const & cdr);

/// @brief Calculate the distance between the carbon and nitrogen of each residue before and after the CDR respectively
core::Real
cdr_CN_anchor_distance(core::pose::Pose const & pose, AntibodyInfoCOP ab_info, CDRNameEnum const & cdr);


}
}


#endif

