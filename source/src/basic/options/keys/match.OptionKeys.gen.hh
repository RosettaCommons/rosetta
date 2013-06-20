// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/match.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_match_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_match_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace match { extern BooleanOptionKey const match; }
namespace match { extern StringOptionKey const lig_name; }
namespace match { extern RealOptionKey const bump_tolerance; }
namespace match { extern FileOptionKey const active_site_definition_by_residue; }
namespace match { extern FileOptionKey const active_site_definition_by_gridlig; }
namespace match { extern FileOptionKey const required_active_site_atom_names; }
namespace match { extern FileOptionKey const grid_boundary; }
namespace match { extern FileOptionKey const geometric_constraint_file; }
namespace match { extern FileOptionKey const scaffold_active_site_residues; }
namespace match { extern FileOptionKey const scaffold_active_site_residues_for_geomcsts; }
namespace match { extern RealOptionKey const euclid_bin_size; }
namespace match { extern RealOptionKey const euler_bin_size; }
namespace match { extern BooleanOptionKey const consolidate_matches; }
namespace match { extern IntegerOptionKey const output_matches_per_group; }
namespace match { extern StringVectorOptionKey const orientation_atoms; }
namespace match { extern StringOptionKey const output_format; }
namespace match { extern StringOptionKey const match_grouper; }
namespace match { extern RealOptionKey const grouper_downstream_rmsd; }
namespace match { extern BooleanOptionKey const output_matchres_only; }
namespace match { extern IntegerVectorOptionKey const geom_csts_downstream_output; }
namespace match { extern BooleanOptionKey const filter_colliding_upstream_residues; }
namespace match { extern RealOptionKey const upstream_residue_collision_tolerance; }
namespace match { extern RealOptionKey const upstream_residue_collision_score_cutoff; }
namespace match { extern RealOptionKey const upstream_residue_collision_Wfa_atr; }
namespace match { extern RealOptionKey const upstream_residue_collision_Wfa_rep; }
namespace match { extern RealOptionKey const upstream_residue_collision_Wfa_sol; }
namespace match { extern BooleanOptionKey const filter_upstream_downstream_collisions; }
namespace match { extern RealOptionKey const updown_collision_tolerance; }
namespace match { extern RealOptionKey const updown_residue_collision_score_cutoff; }
namespace match { extern RealOptionKey const updown_residue_collision_Wfa_atr; }
namespace match { extern RealOptionKey const updown_residue_collision_Wfa_rep; }
namespace match { extern RealOptionKey const updown_residue_collision_Wfa_sol; }
namespace match { extern BooleanOptionKey const define_match_by_single_downstream_positioning; }
namespace match { extern IntegerOptionKey const ligand_rotamer_index; }
namespace match { extern BooleanOptionKey const enumerate_ligand_rotamers; }
namespace match { extern BooleanOptionKey const only_enumerate_non_match_redundant_ligand_rotamers; }
namespace match { extern BooleanOptionKey const dynamic_grid_refinement; }
namespace match { extern BooleanOptionKey const build_round1_hits_twice; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
