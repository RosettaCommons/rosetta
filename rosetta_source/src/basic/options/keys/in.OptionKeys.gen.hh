// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/in.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_in_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_in_OptionKeys_gen_HH

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

namespace in { extern BooleanOptionKey const in; }
namespace in { extern BooleanOptionKey const termini; }
namespace in { extern BooleanOptionKey const ignore_unrecognized_res; }
namespace in { extern BooleanOptionKey const ignore_waters; }
namespace in { extern BooleanOptionKey const add_orbitals; }
namespace in { extern BooleanOptionKey const remember_unrecognized_res; }
namespace in { extern BooleanOptionKey const remember_unrecognized_water; }
namespace in { extern BooleanOptionKey const detect_disulf; }
namespace in { extern FileOptionKey const fix_disulf; }
namespace in { extern BooleanOptionKey const use_stupid_foldtree_format; }
namespace in { extern IntegerVectorOptionKey const target_residues; }
namespace in { extern IntegerVectorOptionKey const replonly_residues; }
namespace in { extern BooleanOptionKey const replonly_loops; }
namespace in { extern BooleanOptionKey const use_database; }
namespace in { extern StringVectorOptionKey const select_structures_from_database; }
namespace in { namespace path { extern PathVectorOptionKey const path; } }
namespace in { namespace path { extern PathVectorOptionKey const fragments; } }
namespace in { namespace path { extern PathVectorOptionKey const pdb; } }
namespace in { namespace path { extern PathVectorOptionKey const database; } }
namespace in { namespace file { extern BooleanOptionKey const file; } }
namespace in { namespace file { extern FileVectorOptionKey const s; } }
namespace in { namespace file { extern FileVectorOptionKey const l; } }
namespace in { namespace file { extern FileVectorOptionKey const list; } }
namespace in { namespace file { extern FileOptionKey const native; } }
namespace in { namespace file { extern FileOptionKey const torsion_bin_probs; } }
namespace in { namespace file { extern FileOptionKey const PCS_frag_cst; } }
namespace in { namespace file { extern FileOptionKey const talos_phi_psi; } }
namespace in { namespace file { extern FileOptionKey const talos_cs; } }
namespace in { namespace file { extern FileOptionKey const ambig_talos_cs_A; } }
namespace in { namespace file { extern FileOptionKey const ambig_talos_cs_B; } }
namespace in { namespace file { extern IntegerVectorOptionKey const native_exclude_res; } }
namespace in { namespace file { extern StringVectorOptionKey const tags; } }
namespace in { namespace file { extern StringVectorOptionKey const user_tags; } }
namespace in { namespace file { extern FileOptionKey const tagfile; } }
namespace in { namespace file { extern FileVectorOptionKey const frag_files; } }
namespace in { namespace file { extern IntegerVectorOptionKey const frag_sizes; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_fa; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_mol; } }
namespace in { namespace file { extern FileVectorOptionKey const extra_res_cen; } }
namespace in { namespace file { extern PathVectorOptionKey const extra_res_path; } }
namespace in { namespace file { extern StringOptionKey const frag3; } }
namespace in { namespace file { extern StringOptionKey const frag9; } }
namespace in { namespace file { extern StringOptionKey const fragA; } }
namespace in { namespace file { extern StringOptionKey const fragB; } }
namespace in { namespace file { extern StringOptionKey const xyz; } }
namespace in { namespace file { extern IntegerOptionKey const fragA_size; } }
namespace in { namespace file { extern BooleanOptionKey const keep_input_scores; } }
namespace in { namespace file { extern BooleanOptionKey const lazy_silent; } }
namespace in { namespace file { extern FileVectorOptionKey const silent; } }
namespace in { namespace file { extern FileVectorOptionKey const atom_tree_diff; } }
namespace in { namespace file { extern StringOptionKey const zip; } }
namespace in { namespace file { extern FileVectorOptionKey const boinc_wu_zip; } }
namespace in { namespace file { extern BooleanOptionKey const fullatom; } }
namespace in { namespace file { extern BooleanOptionKey const centroid_input; } }
namespace in { namespace file { extern BooleanOptionKey const centroid; } }
namespace in { namespace file { extern StringOptionKey const residue_type_set; } }
namespace in { namespace file { extern FileOptionKey const pca; } }
namespace in { namespace file { extern RealOptionKey const silent_energy_cut; } }
namespace in { namespace file { extern FileVectorOptionKey const silent_list; } }
namespace in { namespace file { extern BooleanOptionKey const silent_renumber; } }
namespace in { namespace file { extern BooleanOptionKey const silent_optH; } }
namespace in { namespace file { extern StringOptionKey const silent_struct_type; } }
namespace in { namespace file { extern BooleanOptionKey const silent_read_through_errors; } }
namespace in { namespace file { extern StringOptionKey const silent_score_prefix; } }
namespace in { namespace file { extern IntegerOptionKey const silent_select_random; } }
namespace in { namespace file { extern StringVectorOptionKey const silent_scores_wanted; } }
namespace in { namespace file { extern FileVectorOptionKey const fasta; } }
namespace in { namespace file { extern FileVectorOptionKey const pssm; } }
namespace in { namespace file { extern StringVectorOptionKey const seq; } }
namespace in { namespace file { extern FileOptionKey const checkpoint; } }
namespace in { namespace file { extern FileVectorOptionKey const alignment; } }
namespace in { namespace file { extern FileVectorOptionKey const alignment2; } }
namespace in { namespace file { extern FileOptionKey const rama2b_map; } }
namespace in { namespace file { extern FileOptionKey const psipred_ss2; } }
namespace in { namespace file { extern FileOptionKey const dssp; } }
namespace in { namespace file { extern BooleanOptionKey const fail_on_bad_hbond; } }
namespace in { namespace file { extern FileOptionKey const movemap; } }
namespace in { namespace file { extern BooleanOptionKey const repair_sidechains; } }
namespace in { namespace file { extern BooleanOptionKey const no_binary_dunlib; } }
namespace in { namespace file { extern IntegerOptionKey const extended_pose; } }
namespace in { namespace file { extern FileVectorOptionKey const template_pdb; } }
namespace in { namespace file { extern FileOptionKey const template_silent; } }
namespace in { namespace file { extern FileVectorOptionKey const rdc; } }
namespace in { namespace file { extern FileVectorOptionKey const burial; } }
namespace in { namespace file { extern FileVectorOptionKey const vall; } }
namespace in { namespace file { extern BooleanOptionKey const rescore; } }
namespace in { namespace file { extern StringOptionKey const spanfile; } }
namespace in { namespace file { extern StringOptionKey const lipofile; } }
namespace in { namespace file { extern FileOptionKey const sucker_params; } }
namespace in { namespace file { extern FileOptionKey const fold_tree; } }
namespace in { namespace file { extern BooleanOptionKey const obey_ENDMDL; } }
namespace in { namespace file { extern BooleanOptionKey const new_chain_order; } }
namespace in { namespace file { extern FileOptionKey const ddg_predictions_file; } }
namespace in { namespace rdf { extern BooleanOptionKey const rdf; } }
namespace in { namespace rdf { extern BooleanOptionKey const sep_bb_ss; } }
namespace in { namespace matdes_dock { extern BooleanOptionKey const matdes_dock; } }
namespace in { namespace matdes_dock { extern RealVectorOptionKey const radial_disp; } }
namespace in { namespace matdes_dock { extern RealOptionKey const neg_r; } }
namespace in { namespace matdes_dock { extern RealVectorOptionKey const angle; } }
namespace in { namespace matdes_dock { extern BooleanOptionKey const dump_pdb; } }
namespace in { namespace matdes_dock { extern BooleanOptionKey const dump_chainA_only; } }
namespace in { namespace matdes_dock { extern StringOptionKey const pdbID; } }
namespace in { namespace matdes_dock { extern StringOptionKey const prefix; } }
namespace in { namespace matdes_dock { extern IntegerOptionKey const num_subs_building_block; } }
namespace in { namespace matdes_dock { extern IntegerOptionKey const num_subs_total; } }
namespace in { namespace matdes_design { extern BooleanOptionKey const matdes_design; } }
namespace in { namespace matdes_design { extern RealOptionKey const contact_dist; } }
namespace in { namespace matdes_design { extern RealOptionKey const grid_size_angle; } }
namespace in { namespace matdes_design { extern RealOptionKey const grid_size_radius; } }
namespace in { namespace matdes_design { extern IntegerOptionKey const grid_nsamp_angle; } }
namespace in { namespace matdes_design { extern IntegerOptionKey const grid_nsamp_radius; } }
namespace in { namespace matdes_design { extern IntegerOptionKey const num_subs_building_block; } }
namespace in { namespace matdes_design { extern RealOptionKey const fav_nat_bonus; } }
namespace in { namespace matdes_design { extern IntegerVectorOptionKey const revert_pos; } }
namespace in { namespace matdes_design { extern StringVectorOptionKey const revert_ids; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
