// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/relax.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_relax_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_relax_OptionKeys_gen_HH

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

namespace relax { extern BooleanOptionKey const relax; }
namespace relax { extern BooleanOptionKey const fast; }
namespace relax { extern BooleanOptionKey const thorough; }
namespace relax { extern BooleanOptionKey const membrane; }
namespace relax { extern IntegerOptionKey const default_repeats; }
namespace relax { extern BooleanOptionKey const ramady; }
namespace relax { extern RealOptionKey const ramady_rms_limit; }
namespace relax { extern RealOptionKey const ramady_cutoff; }
namespace relax { extern IntegerOptionKey const ramady_max_rebuild; }
namespace relax { extern BooleanOptionKey const ramady_force; }
namespace relax { extern FileOptionKey const script; }
namespace relax { extern IntegerOptionKey const script_max_accept; }
namespace relax { extern BooleanOptionKey const superimpose_to_native; }
namespace relax { extern FileOptionKey const superimpose_to_file; }
namespace relax { extern BooleanOptionKey const constrain_relax_to_native_coords; }
namespace relax { extern BooleanOptionKey const constrain_relax_to_start_coords; }
namespace relax { extern RealOptionKey const sc_cst_maxdist; }
namespace relax { extern BooleanOptionKey const limit_aroma_chi2; }
namespace relax { extern BooleanOptionKey const bb_move; }
namespace relax { extern BooleanOptionKey const chi_move; }
namespace relax { extern BooleanOptionKey const jump_move; }
namespace relax { extern BooleanOptionKey const minimize_bond_lengths; }
namespace relax { extern BooleanOptionKey const minimize_bond_angles; }
namespace relax { extern BooleanOptionKey const minimize_mainchain_bond_lengths; }
namespace relax { extern BooleanOptionKey const minimize_mainchain_bond_angles; }
namespace relax { extern StringOptionKey const min_type; }
namespace relax { extern BooleanOptionKey const cartesian; }
namespace relax { extern RealOptionKey const chainbreak_weight; }
namespace relax { extern RealOptionKey const linear_chainbreak_weight; }
namespace relax { extern RealOptionKey const overlap_chainbreak_weight; }
namespace relax { extern BooleanOptionKey const classic; }
namespace relax { extern FileOptionKey const sequence_file; }
namespace relax { extern BooleanOptionKey const quick; }
namespace relax { extern BooleanOptionKey const sequence; }
namespace relax { extern IntegerOptionKey const minirelax_repeats; }
namespace relax { extern RealOptionKey const minirelax_sdev; }
namespace relax { extern BooleanOptionKey const wobblemoves; }
namespace relax { extern FileOptionKey const constrain_relax_segments; }
namespace relax { extern RealOptionKey const coord_cst_width; }
namespace relax { extern RealOptionKey const coord_cst_stdev; }
namespace relax { extern BooleanOptionKey const ramp_constraints; }
namespace relax { extern RealOptionKey const energycut; }
namespace relax { extern BooleanOptionKey const mini; }
namespace relax { extern IntegerOptionKey const stage1_ramp_cycles; }
namespace relax { extern IntegerOptionKey const stage1_ramp_inner_cycles; }
namespace relax { extern IntegerOptionKey const stage2_repack_period; }
namespace relax { extern IntegerOptionKey const stage2_cycles; }
namespace relax { extern RealOptionKey const min_tolerance; }
namespace relax { extern IntegerOptionKey const stage3_cycles; }
namespace relax { extern RealOptionKey const cycle_ratio; }
namespace relax { extern RealOptionKey const filter_stage2_beginning; }
namespace relax { extern RealOptionKey const filter_stage2_quarter; }
namespace relax { extern RealOptionKey const filter_stage2_half; }
namespace relax { extern RealOptionKey const filter_stage2_end; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
