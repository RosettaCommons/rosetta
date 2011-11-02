// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/evaluation.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_evaluation_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_evaluation_OptionKeys_gen_HH

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

namespace evaluation { extern BooleanOptionKey const evaluation; }
namespace evaluation { extern FileVectorOptionKey const rmsd_target; }
namespace evaluation { extern StringVectorOptionKey const rmsd_column; }
namespace evaluation { extern FileVectorOptionKey const rmsd_select; }
namespace evaluation { extern FileVectorOptionKey const align_rmsd_target; }
namespace evaluation { extern FileVectorOptionKey const structural_similarity; }
namespace evaluation { extern BooleanOptionKey const contact_map; }
namespace evaluation { extern StringVectorOptionKey const jscore_evaluator; }
namespace evaluation { extern StringVectorOptionKey const align_rmsd_column; }
namespace evaluation { extern FileVectorOptionKey const align_rmsd_fns; }
namespace evaluation { extern StringOptionKey const align_rmsd_format; }
namespace evaluation { extern StringOptionKey const predicted_burial_fn; }
namespace evaluation { extern FileOptionKey const pool; }
namespace evaluation { extern FileVectorOptionKey const rmsd; }
namespace evaluation { extern BooleanOptionKey const gdtmm; }
namespace evaluation { extern BooleanOptionKey const score_with_rmsd; }
namespace evaluation { extern FileVectorOptionKey const constraints; }
namespace evaluation { extern FileVectorOptionKey const constraints_column; }
namespace evaluation { extern FileVectorOptionKey const combined_constraints; }
namespace evaluation { extern FileVectorOptionKey const combined_constraints_column; }
namespace evaluation { extern IntegerOptionKey const combine_statistics; }
namespace evaluation { extern StringVectorOptionKey const chemical_shifts; }
namespace evaluation { extern StringOptionKey const sparta_dir; }
namespace evaluation { extern StringVectorOptionKey const cam_shifts; }
namespace evaluation { extern StringVectorOptionKey const pales; }
namespace evaluation { extern FileVectorOptionKey const extra_score; }
namespace evaluation { extern FileVectorOptionKey const extra_score_patch; }
namespace evaluation { extern StringVectorOptionKey const extra_score_column; }
namespace evaluation { extern FileVectorOptionKey const extra_score_select; }
namespace evaluation { extern FileVectorOptionKey const rdc_select; }
namespace evaluation { extern FileVectorOptionKey const rdc_target; }
namespace evaluation { extern BooleanOptionKey const symmetric_rmsd; }
namespace evaluation { extern StringVectorOptionKey const rdc_column; }
namespace evaluation { extern StringVectorOptionKey const rdc; }
namespace evaluation { extern BooleanOptionKey const jump_nr; }
namespace evaluation { extern IntegerVectorOptionKey const score_exclude_res; }
namespace evaluation { extern IntegerOptionKey const score_sscore_short_helix; }
namespace evaluation { extern IntegerOptionKey const score_sscore_maxloop; }
namespace evaluation { extern BooleanOptionKey const rpf; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
