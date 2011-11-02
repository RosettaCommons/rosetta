// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/willmatch.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_willmatch_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_willmatch_OptionKeys_gen_HH

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

namespace willmatch { extern BooleanOptionKey const willmatch; }
namespace willmatch { extern RealOptionKey const arg_dun_th; }
namespace willmatch { extern RealOptionKey const asp_dun_th; }
namespace willmatch { extern RealOptionKey const glu_dun_th; }
namespace willmatch { extern RealOptionKey const lys_dun_th; }
namespace willmatch { extern BooleanOptionKey const usecache; }
namespace willmatch { extern StringVectorOptionKey const write_reduced_matchset; }
namespace willmatch { extern RealOptionKey const interface_size; }
namespace willmatch { extern RealOptionKey const max_dis_any; }
namespace willmatch { extern RealOptionKey const max_dis_all; }
namespace willmatch { extern RealOptionKey const max_dis_hb; }
namespace willmatch { extern RealOptionKey const min_dis_hb; }
namespace willmatch { extern RealOptionKey const max_dis_hb_colinear; }
namespace willmatch { extern RealOptionKey const max_dis_metal; }
namespace willmatch { extern RealOptionKey const max_ang_metal; }
namespace willmatch { extern RealOptionKey const clash_dis; }
namespace willmatch { extern RealOptionKey const c2_linker_dist; }
namespace willmatch { extern RealOptionKey const identical_match_dis; }
namespace willmatch { extern RealOptionKey const chi1_increment; }
namespace willmatch { extern RealOptionKey const chi2_increment; }
namespace willmatch { extern RealOptionKey const c2_symm_increment; }
namespace willmatch { extern RealOptionKey const cb_sasa_thresh; }
namespace willmatch { extern BooleanOptionKey const design_interface; }
namespace willmatch { extern FileOptionKey const chilist; }
namespace willmatch { extern FileOptionKey const fixed_res; }
namespace willmatch { extern FileOptionKey const native1; }
namespace willmatch { extern FileOptionKey const native2; }
namespace willmatch { extern FileOptionKey const exclude_res1; }
namespace willmatch { extern FileOptionKey const exclude_res2; }
namespace willmatch { extern FileOptionKey const taglist; }
namespace willmatch { extern IntegerVectorOptionKey const residues; }
namespace willmatch { extern BooleanOptionKey const symmetry_d2; }
namespace willmatch { extern BooleanOptionKey const symmetry_c2_dock; }
namespace willmatch { extern IntegerVectorOptionKey const splitwork; }
namespace willmatch { extern BooleanOptionKey const exclude_ala; }
namespace willmatch { extern RealOptionKey const match_overlap_dis; }
namespace willmatch { extern RealOptionKey const match_overlap_ang; }
namespace willmatch { extern IntegerVectorOptionKey const forbid_residues; }
namespace willmatch { extern RealVectorOptionKey const poi; }
namespace willmatch { extern RealOptionKey const poidis; }
namespace willmatch { extern BooleanOptionKey const homodimer; }
namespace willmatch { extern RealOptionKey const fa_dun_thresh; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
