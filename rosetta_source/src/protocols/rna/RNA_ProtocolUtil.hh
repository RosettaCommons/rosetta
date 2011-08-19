// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_ProtocolUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_ProtocolUtil_hh
#define INCLUDED_protocols_rna_RNA_ProtocolUtil_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/io/ozstream.fwd.hh>

// Utility headers
// Auto-header: duplicate removed #include <utility/io/ozstream.fwd.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <vector>

//Auto Headers
#include <ObjexxFCL/FArray1D.fwd.hh>



namespace protocols {
namespace rna {

typedef  numeric::xyzVector< core::Real >  Vector;

void
figure_out_reasonable_rna_fold_tree( core::pose::Pose & pose );

void
get_base_pairing_info( core::pose::Pose const & pose,
											 core::Size const & seqpos,
											 char & secstruct,
											 ObjexxFCL::FArray1D <bool> & edge_is_base_pairing );

void
create_rna_vall_torsions( core::pose::Pose & pose, utility::io::ozstream & torsions_out);

void
create_rna_vall_torsions( core::pose::Pose & pose, std::string const outfile );

core::Real
get_o1p_o2p_sign( core::pose::Pose const & pose );

void
ensure_phosphate_nomenclature_matches_mini( core::pose::Pose & pose );

core::Real
get_o1p_o2p_sign_parin( core::pose::Pose const & pose , core::Size res_num);

void
assert_phosphate_nomenclature_matches_mini( core::pose::Pose const & pose);

void
ensure_phosphate_nomenclature_matches_mini_parin ( core::pose::Pose & pose);

void
export_packer_results(
	utility::vector1< std::pair< core::Real, std::string > > & results,
	utility::vector1< core::pose::PoseOP > pose_list,
	core::scoring::ScoreFunctionOP & scorefxn,
	std::string const & outfile,
	bool const dumo = false );

void
check_base_pair( core::pose::Pose & pose, ObjexxFCL::FArray1D_int & struct_type );

void
setup_base_pair_constraints(
														core::pose::Pose & pose,
														utility::vector1< std::pair< core::Size, core::Size > > const &  pairings );

} //rna
} // protocols

#endif
