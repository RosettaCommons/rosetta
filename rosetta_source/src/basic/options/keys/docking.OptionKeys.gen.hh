// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/docking.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_docking_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_docking_OptionKeys_gen_HH

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

namespace docking { extern BooleanOptionKey const kick_relax; }
namespace docking { extern BooleanOptionKey const docking; }
namespace docking { extern BooleanOptionKey const view; }
namespace docking { extern BooleanOptionKey const no_filters; }
namespace docking { extern StringVectorOptionKey const design_chains; }
namespace docking { extern FileOptionKey const recover_sidechains; }
namespace docking { extern StringOptionKey const partners; }
namespace docking { extern BooleanOptionKey const docking_local_refine; }
namespace docking { extern BooleanOptionKey const low_res_protocol_only; }
namespace docking { extern BooleanOptionKey const randomize1; }
namespace docking { extern BooleanOptionKey const randomize2; }
namespace docking { extern BooleanOptionKey const spin; }
namespace docking { extern RealVectorOptionKey const dock_pert; }
namespace docking { extern RealOptionKey const uniform_trans; }
namespace docking { extern BooleanOptionKey const center_at_interface; }
namespace docking { extern IntegerOptionKey const docking_centroid_outer_cycles; }
namespace docking { extern IntegerOptionKey const docking_centroid_inner_cycles; }
namespace docking { extern BooleanOptionKey const dock_min; }
namespace docking { extern StringOptionKey const flexible_bb_docking; }
namespace docking { extern RealOptionKey const flexible_bb_docking_interface_dist; }
namespace docking { extern StringOptionKey const ensemble1; }
namespace docking { extern StringOptionKey const ensemble2; }
namespace docking { extern RealOptionKey const dock_mcm_trans_magnitude; }
namespace docking { extern RealOptionKey const dock_mcm_rot_magnitude; }
namespace docking { extern BooleanOptionKey const dock_rtmin; }
namespace docking { extern BooleanOptionKey const sc_min; }
namespace docking { extern BooleanOptionKey const norepack1; }
namespace docking { extern BooleanOptionKey const norepack2; }
namespace docking { extern BooleanOptionKey const dock_ppk; }
namespace docking { extern IntegerOptionKey const max_repeats; }
namespace docking { extern RealVectorOptionKey const dock_lowres_filter; }
namespace docking { extern IntegerVectorOptionKey const multibody; }
namespace docking { extern BooleanOptionKey const ignore_default_docking_task; }
namespace docking { extern StringOptionKey const low_patch; }
namespace docking { extern StringOptionKey const high_patch; }
namespace docking { extern StringOptionKey const high_min_patch; }
namespace docking { extern StringOptionKey const pack_patch; }
namespace docking { extern BooleanOptionKey const use_legacy_protocol; }
namespace docking { namespace ligand { extern BooleanOptionKey const ligand; } }
namespace docking { namespace ligand { extern StringOptionKey const protocol; } }
namespace docking { namespace ligand { extern BooleanOptionKey const soft_rep; } }
namespace docking { namespace ligand { extern BooleanOptionKey const tweak_sxfn; } }
namespace docking { namespace ligand { extern BooleanOptionKey const old_estat; } }
namespace docking { namespace ligand { extern BooleanOptionKey const random_conformer; } }
namespace docking { namespace ligand { extern IntegerOptionKey const improve_orientation; } }
namespace docking { namespace ligand { extern BooleanOptionKey const mutate_same_name3; } }
namespace docking { namespace ligand { extern RealOptionKey const min_rms; } }
namespace docking { namespace ligand { extern IntegerOptionKey const max_poses; } }
namespace docking { namespace ligand { extern BooleanOptionKey const minimize_ligand; } }
namespace docking { namespace ligand { extern RealOptionKey const harmonic_torsions; } }
namespace docking { namespace ligand { extern BooleanOptionKey const use_ambig_constraints; } }
namespace docking { namespace ligand { extern IntegerOptionKey const shear_moves; } }
namespace docking { namespace ligand { extern BooleanOptionKey const minimize_backbone; } }
namespace docking { namespace ligand { extern RealOptionKey const harmonic_Calphas; } }
namespace docking { namespace ligand { extern RealOptionKey const tether_ligand; } }
namespace docking { namespace ligand { extern RealVectorOptionKey const start_from; } }
namespace docking { namespace ligand { extern StringOptionKey const option_file; } }
namespace docking { namespace ligand { extern BooleanOptionKey const rescore; } }
namespace docking { namespace ligand { namespace grid { extern BooleanOptionKey const grid; } } }
namespace docking { namespace ligand { namespace grid { extern FileOptionKey const grid_kin; } } }
namespace docking { namespace ligand { namespace grid { extern FileOptionKey const grid_map; } } }
namespace docking { namespace symmetry { extern BooleanOptionKey const symmetry; } }
namespace docking { namespace symmetry { extern BooleanOptionKey const minimize_backbone; } }
namespace docking { namespace symmetry { extern BooleanOptionKey const minimize_sidechains; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
