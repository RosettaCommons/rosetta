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

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>


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

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
