// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_denovo_util_HH
#define INCLUDED_protocols_rna_denovo_util_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/pose/rna/BasePairStep.fwd.hh>
#include <core/pose/rna/RNA_BaseDoubletClasses.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <map>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


namespace protocols {
namespace rna {
namespace denovo {

typedef  numeric::xyzVector< core::Length >  Vector;

void
create_rna_vall_torsions( core::pose::Pose & pose,
	utility::io::ozstream & torsions_out,
	utility::vector1 <core::Size> const & exclude_res_list );

void
create_rna_vall_torsions( core::pose::Pose & pose,
	std::string const & outfile,
	utility::vector1 <core::Size> const & exclude_res_list );

void
ensure_phosphate_nomenclature_matches_mini( core::pose::Pose & pose );

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
setup_coarse_chainbreak_constraints( core::pose::Pose & pose, core::Size const & n );
void
print_internal_coords( core::pose::Pose const & pose );

bool
possible_root( core::kinematics::FoldTree const & f, core::Size const & n );

inline bool is_num_in_list ( core::Size const i,
	utility::vector1 <core::Size> const & list )
{
	return std::find(list.begin(), list.end(), i)!=list.end();
}

bool
let_rigid_body_jumps_move( core::kinematics::MoveMap & movemap,
	core::pose::Pose const & pose,
	bool const & move_first_rigid_body  = false,
	bool const & allow_first_rigid_body = false );

void
translate_virtual_anchor_to_first_rigid_body( core::pose::Pose & pose );

bool
involved_in_phosphate_torsion( std::string const & atomname );

void
figure_out_base_pair_partner( core::pose::Pose & pose, std::map< core::Size, core::Size > & partner,
	bool const strict = true );

utility::vector1< core::pose::rna::BasePair >
classify_base_pairs_lores( core::pose::Pose const & pose );

void
print_hbonds( core::pose::Pose & pose );

core::Size
virtualize_bulges( core::pose::Pose & input_pose,
	utility::vector1< core::Size > const & in_allow_bulge_res_list,
	core::scoring::ScoreFunctionCOP const & scorefxn,
	std::string const & tag,
	bool const allow_pre_virtualize,
	bool const allow_consecutive_bulges,
	bool const verbose );

core::scoring::ScoreFunctionOP
get_rna_hires_scorefxn( bool const & include_protein_terms=false );

utility::vector1< core::Size >
get_moving_res( core::pose::Pose const & pose,
	core::pose::toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

void delete_non_protein_from_pose( core::pose::Pose & pose );

utility::vector1< core::Size >
get_residues_within_dist_of_RNA(
	core::pose::Pose const & pose,
	core::Real const & dist_cutoff );


core::kinematics::FoldTree
get_rnp_docking_fold_tree( core::pose::Pose const & pose,
	bool const & with_density = false );


} //denovo
} //rna
} //protocols

#endif
