// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_util_HH
#define INCLUDED_protocols_magnesium_util_HH

#include <protocols/magnesium/params.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <numeric/UniformRotationSampler.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <map>

namespace protocols {
namespace magnesium {

void
fixup_magnesiums( core::pose::Pose & pose );

void
hydrate_magnesiums( core::pose::Pose & pose,
	bool use_virtual_waters_as_placeholders = true,
	bool test_all_mg_hydration_frames = false );

core::scoring::ScoreFunctionOP
get_mg_scorefxn();

numeric::UniformRotationSamplerCOP
get_water_uniform_rotation_sampler();

numeric::UniformRotationSamplerCOP
get_octahedral_uniform_rotation_sampler( core::Real const rotstep = 2.5,
	bool const remove_redundant = false );

core::conformation::ResidueOP
get_useful_HOH_coords( core::Vector & Oc, core::Vector & OH1c, core::Vector & OH2c,
	core::chemical::ResidueTypeSet const & residue_set );

void
add_single_magnesium( core::pose::Pose & pose );

void
strip_out_magnesiums( core::pose::Pose & pose );

std::map< core::Size, utility::vector1< core::Size > >
define_mg_water_map( core::pose::Pose const & pose );

utility::vector1< std::pair< core::Size, core::Size > >
get_mg_water_pairs( core::pose::Pose const & pose,
	bool const exclude_virtual_waters = true ) ;

utility::vector1< std::pair< core::Size, core::Size > >
get_mg_water_pairs( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & mg_res,
	bool const exclude_virtual_waters = true ) ;

utility::vector1< core::id::AtomID >
get_mg_ligands( core::pose::Pose const & pose, core::Size const i,
	bool const filter_for_acceptors = true,
	bool const exclude_virtual_waters = true );

utility::vector1< core::id::AtomID >
filter_acceptor_ligands( core::pose::Pose const & pose, utility::vector1< core::id::AtomID > const & ligands );

core::Size
instantiate_water_at_octahedral_vertex( core::pose::Pose & pose,
	core::Size const mg_res,
	core::Size const n /* 1 ... 6*/,
	core::Distance const hoh_distance = MG_HOH_DISTANCE,
	bool const replace_residue  = false,
	bool const virtual_water  = false );

void
update_jump_atoms_for_mg_bound_water( core::pose::Pose & pose, core::Size const n );

core::Size
append_mg_bound_water(  core::pose::Pose & pose,
	core::conformation::Residue const & rsd,
	core::Size const mg_res );

utility::vector1< core::Size >
find_bound_waters_that_are_daughters_in_fold_tree( core::pose::Pose const & pose, core::Size const mg_res );

core::conformation::ResidueOP
get_mg_rsd();

utility::vector1< core::Size >
get_mg_res( core::pose::Pose const & pose );

utility::vector1< core::Size >
get_water_res( core::pose::Pose const & pose );

void
remove_waters_except_mg_bound( core::pose::Pose & pose,
	utility::vector1< std::pair< core::Size, core::Size > > const & mg_water_pairs );

void
remove_mg_bound_waters( core::pose::Pose & pose, utility::vector1< core::Size > const & mg_res, bool const leave_other_waters = false );


void
update_numbers_in_pdb_info( core::pose::Pose & pose, bool const reset_waters = false );

utility::vector1< core::Size >
pdbslice( core::pose::Pose & pose, core::Size const center_res, core::Distance distance_cutoff = 12.0 );

void
get_hydration_stats( core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	utility::vector1< core::Size > const & pdb_mg_res_list_in,
	std::string const & outfile );

} //magnesium
} //protocols



#endif
