// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/ddg.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_ddg_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_ddg_OptionKeys_gen_HH

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

namespace ddg { extern BooleanOptionKey const ddg; }
namespace ddg { extern BooleanOptionKey const avg_rot_cst_enrg; }
namespace ddg { extern BooleanOptionKey const use_bound_cst; }
namespace ddg { extern RealOptionKey const cap_rot_cst_enrg; }
namespace ddg { extern BooleanOptionKey const opt_input_structure; }
namespace ddg { extern BooleanOptionKey const pack_until_converge; }
namespace ddg { extern BooleanOptionKey const no_constraints; }
namespace ddg { extern RealOptionKey const apply_rot_cst_to_mutation_region_only; }
namespace ddg { extern RealOptionKey const rot_cst_enrg_cutoff; }
namespace ddg { extern BooleanOptionKey const use_rotamer_constraints_to_native; }
namespace ddg { extern BooleanOptionKey const single_res_scoring; }
namespace ddg { extern BooleanOptionKey const downweight_by_sasa; }
namespace ddg { extern BooleanOptionKey const global; }
namespace ddg { extern BooleanOptionKey const exclude_solvent_exposed_res; }
namespace ddg { extern RealOptionKey const radius; }
namespace ddg { extern StringOptionKey const wt; }
namespace ddg { extern StringOptionKey const mut; }
namespace ddg { extern BooleanOptionKey const suppress_checkpointing; }
namespace ddg { extern BooleanOptionKey const wt_only; }
namespace ddg { extern BooleanOptionKey const mut_only; }
namespace ddg { extern BooleanOptionKey const output_silent; }
namespace ddg { extern StringOptionKey const minimization_scorefunction; }
namespace ddg { extern StringOptionKey const minimization_patch; }
namespace ddg { extern BooleanOptionKey const min_cst; }
namespace ddg { extern IntegerOptionKey const lowest_x_decoys; }
namespace ddg { extern BooleanOptionKey const local_opt_only; }
namespace ddg { extern BooleanOptionKey const print_per_res_diff; }
namespace ddg { extern BooleanOptionKey const mean; }
namespace ddg { extern BooleanOptionKey const min; }
namespace ddg { extern BooleanOptionKey const rb_restrict_to_mutsite_nbrs; }
namespace ddg { extern BooleanOptionKey const no_bb_movement; }
namespace ddg { extern BooleanOptionKey const initial_repack; }
namespace ddg { extern StringOptionKey const rb_file; }
namespace ddg { extern IntegerOptionKey const interface_ddg; }
namespace ddg { extern RealOptionKey const ens_variation; }
namespace ddg { extern BooleanOptionKey const sc_min_only; }
namespace ddg { extern StringOptionKey const min_cst_weights; }
namespace ddg { extern RealOptionKey const opt_radius; }
namespace ddg { extern StringOptionKey const output_dir; }
namespace ddg { extern StringOptionKey const last_accepted_pose_dir; }
namespace ddg { extern BooleanOptionKey const min_with_cst; }
namespace ddg { extern RealOptionKey const temperature; }
namespace ddg { extern BooleanOptionKey const ramp_repulsive; }
namespace ddg { extern StringOptionKey const mut_file; }
namespace ddg { extern StringOptionKey const out_pdb_prefix; }
namespace ddg { extern RealOptionKey const constraint_weight; }
namespace ddg { extern RealOptionKey const harmonic_ca_tether; }
namespace ddg { extern IntegerOptionKey const iterations; }
namespace ddg { extern StringOptionKey const out; }
namespace ddg { extern BooleanOptionKey const debug_output; }
namespace ddg { extern BooleanOptionKey const dump_pdbs; }
namespace ddg { extern StringOptionKey const weight_file; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
