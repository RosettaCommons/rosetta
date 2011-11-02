// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/frags.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_frags_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_frags_OptionKeys_gen_HH

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

namespace frags { extern BooleanOptionKey const frags; }
namespace frags { extern BooleanOptionKey const filter_JC; }
namespace frags { extern BooleanOptionKey const bounded_protocol; }
namespace frags { extern BooleanOptionKey const keep_all_protocol; }
namespace frags { extern BooleanOptionKey const quota_protocol; }
namespace frags { extern BooleanOptionKey const nonlocal_pairs_protocol; }
namespace frags { extern BooleanOptionKey const p_value_selection; }
namespace frags { extern IntegerOptionKey const n_frags; }
namespace frags { extern FileOptionKey const allowed_pdb; }
namespace frags { extern StringVectorOptionKey const ss_pred; }
namespace frags { extern FileOptionKey const denied_pdb; }
namespace frags { extern IntegerVectorOptionKey const frag_sizes; }
namespace frags { extern BooleanOptionKey const write_ca_coordinates; }
namespace frags { extern BooleanOptionKey const annotate; }
namespace frags { extern IntegerOptionKey const nr_large_copies; }
namespace frags { extern IntegerOptionKey const n_candidates; }
namespace frags { extern BooleanOptionKey const write_rama_tables; }
namespace frags { extern RealOptionKey const rama_C; }
namespace frags { extern RealOptionKey const rama_B; }
namespace frags { extern RealOptionKey const sigmoid_cs_A; }
namespace frags { extern RealOptionKey const sigmoid_cs_B; }
namespace frags { extern RealOptionKey const seqsim_H; }
namespace frags { extern RealOptionKey const seqsim_E; }
namespace frags { extern RealOptionKey const seqsim_L; }
namespace frags { extern RealOptionKey const rama_norm; }
namespace frags { extern StringOptionKey const describe_fragments; }
namespace frags { namespace scoring { extern BooleanOptionKey const scoring; } }
namespace frags { namespace scoring { extern FileOptionKey const config; } }
namespace frags { namespace scoring { extern StringOptionKey const profile_score; } }
namespace frags { namespace scoring { extern FileVectorOptionKey const predicted_secondary; } }
namespace frags { namespace picking { extern BooleanOptionKey const picking; } }
namespace frags { namespace picking { extern StringOptionKey const selecting_rule; } }
namespace frags { namespace picking { extern StringOptionKey const selecting_scorefxn; } }
namespace frags { namespace picking { extern FileOptionKey const quota_config_file; } }
namespace frags { namespace picking { extern IntegerVectorOptionKey const query_pos; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const nonlocal; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const relax_input; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const relax_input_with_coordinate_constraints; } }
namespace frags { namespace nonlocal { extern IntegerOptionKey const relax_frags_repeats; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const single_chain; } }
namespace frags { namespace nonlocal { extern IntegerOptionKey const min_seq_sep; } }
namespace frags { namespace nonlocal { extern IntegerOptionKey const ca_dist; } }
namespace frags { namespace nonlocal { extern IntegerOptionKey const min_contacts_per_res; } }
namespace frags { namespace nonlocal { extern RealOptionKey const max_ddg_score; } }
namespace frags { namespace nonlocal { extern RealOptionKey const max_rmsd_after_relax; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const output_frags_pdbs; } }
namespace frags { namespace nonlocal { extern BooleanOptionKey const output_idealized; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
