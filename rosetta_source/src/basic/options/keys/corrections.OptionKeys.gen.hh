// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/corrections.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_corrections_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_corrections_OptionKeys_gen_HH

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

namespace corrections { extern BooleanOptionKey const corrections; }
namespace corrections { extern BooleanOptionKey const correct; }
namespace corrections { namespace score { extern BooleanOptionKey const score; } }
namespace corrections { namespace score { extern BooleanOptionKey const no_his_his_pairE; } }
namespace corrections { namespace score { extern BooleanOptionKey const hbond_His_Phil_fix; } }
namespace corrections { namespace score { extern BooleanOptionKey const helix_hb_06_2009; } }
namespace corrections { namespace score { extern BooleanOptionKey const use_incorrect_hbond_deriv; } }
namespace corrections { namespace score { extern StringOptionKey const p_aa_pp; } }
namespace corrections { namespace score { extern BooleanOptionKey const p_aa_pp_nogridshift; } }
namespace corrections { namespace score { extern BooleanOptionKey const rama_not_squared; } }
namespace corrections { namespace score { extern FileOptionKey const rama_map; } }
namespace corrections { namespace score { extern BooleanOptionKey const dun10; } }
namespace corrections { namespace score { extern StringOptionKey const dun10_dir; } }
namespace corrections { namespace score { extern BooleanOptionKey const dun08; } }
namespace corrections { namespace score { extern StringOptionKey const dun08_dir; } }
namespace corrections { namespace score { extern StringOptionKey const dun02_file; } }
namespace corrections { namespace score { extern StringOptionKey const ch_o_bond_potential; } }
namespace corrections { namespace score { extern BooleanOptionKey const hack_elec_co_only; } }
namespace corrections { namespace score { extern RealOptionKey const lj_hbond_hdis; } }
namespace corrections { namespace score { extern FileOptionKey const PB_potential_file; } }
namespace corrections { namespace score { extern BooleanOptionKey const PB_sidechain_only; } }
namespace corrections { namespace score { extern IntegerVectorOptionKey const PB_score_residue_range; } }
namespace corrections { namespace score { extern IntegerVectorOptionKey const PB_revamp_near_chain; } }
namespace corrections { namespace score { extern RealOptionKey const PB_potential_cap; } }
namespace corrections { namespace score { extern RealOptionKey const lj_hbond_OH_donor_dis; } }
namespace corrections { namespace score { extern BooleanOptionKey const score12prime; } }
namespace corrections { namespace score { extern BooleanOptionKey const hb_sp2_chipen; } }
namespace corrections { namespace score { extern RealOptionKey const hb_sp2_amp; } }
namespace corrections { namespace score { extern RealOptionKey const hb_sp2_peak_heigh_above_trough; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const chemical; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const icoor_05_2009; } }
namespace corrections { namespace chemical { extern BooleanOptionKey const parse_charge; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
