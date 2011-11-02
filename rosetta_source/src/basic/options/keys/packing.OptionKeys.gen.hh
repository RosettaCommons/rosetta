// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/packing.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_packing_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_packing_OptionKeys_gen_HH

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

namespace packing { extern BooleanOptionKey const packing; }
namespace packing { extern BooleanOptionKey const repack_only; }
namespace packing { extern BooleanOptionKey const prevent_repacking; }
namespace packing { extern IntegerOptionKey const ndruns; }
namespace packing { extern BooleanOptionKey const soft_rep_design; }
namespace packing { extern BooleanOptionKey const use_electrostatic_repulsion; }
namespace packing { extern BooleanOptionKey const dump_rotamer_sets; }
namespace packing { extern RealOptionKey const dunbrack_prob_buried; }
namespace packing { extern RealOptionKey const dunbrack_prob_nonburied; }
namespace packing { extern RealOptionKey const dunbrack_prob_nonburied_semirotameric; }
namespace packing { extern BooleanOptionKey const no_optH; }
namespace packing { extern BooleanOptionKey const optH_MCA; }
namespace packing { extern BooleanOptionKey const pack_missing_sidechains; }
namespace packing { extern BooleanOptionKey const preserve_c_beta; }
namespace packing { extern BooleanOptionKey const flip_HNQ; }
namespace packing { extern IntegerVectorOptionKey const fix_his_tautomer; }
namespace packing { extern BooleanOptionKey const print_pymol_selection; }
namespace packing { namespace ex1 { extern BooleanOptionKey const ex1; } }
namespace packing { namespace ex1 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex1 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex2 { extern BooleanOptionKey const ex2; } }
namespace packing { namespace ex2 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex3 { extern BooleanOptionKey const ex3; } }
namespace packing { namespace ex3 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex3 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex4 { extern BooleanOptionKey const ex4; } }
namespace packing { namespace ex4 { extern IntegerOptionKey const level; } }
namespace packing { namespace ex4 { extern BooleanOptionKey const operate; } }
namespace packing { namespace ex1aro { extern BooleanOptionKey const ex1aro; } }
namespace packing { namespace ex1aro { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2aro { extern BooleanOptionKey const ex2aro; } }
namespace packing { namespace ex2aro { extern IntegerOptionKey const level; } }
namespace packing { namespace ex1aro_exposed { extern BooleanOptionKey const ex1aro_exposed; } }
namespace packing { namespace ex1aro_exposed { extern IntegerOptionKey const level; } }
namespace packing { namespace ex2aro_exposed { extern BooleanOptionKey const ex2aro_exposed; } }
namespace packing { namespace ex2aro_exposed { extern IntegerOptionKey const level; } }
namespace packing { namespace exdna { extern BooleanOptionKey const exdna; } }
namespace packing { namespace exdna { extern IntegerOptionKey const level; } }
namespace packing { extern IntegerOptionKey const extrachi_cutoff; }
namespace packing { extern FileVectorOptionKey const resfile; }
namespace packing { extern RealOptionKey const outeriterations_scaling; }
namespace packing { extern RealOptionKey const inneriterations_scaling; }
namespace packing { extern BooleanOptionKey const explicit_h2o; }
namespace packing { extern StringVectorOptionKey const adducts; }
namespace packing { extern BooleanOptionKey const solvate; }
namespace packing { extern BooleanOptionKey const use_input_sc; }
namespace packing { extern FileVectorOptionKey const unboundrot; }
namespace packing { extern RealOptionKey const max_rotbump_energy; }
namespace packing { extern BooleanOptionKey const lazy_ig; }
namespace packing { extern BooleanOptionKey const double_lazy_ig; }
namespace packing { extern IntegerOptionKey const double_lazy_ig_mem_limit; }
namespace packing { extern IntegerOptionKey const linmem_ig; }
namespace packing { extern IntegerOptionKey const multi_cool_annealer; }
namespace packing { extern RealVectorOptionKey const minpack_temp_schedule; }
namespace packing { extern IntegerOptionKey const minpack_inner_iteration_scale; }
namespace packing { extern BooleanOptionKey const minpack_disable_bumpcheck; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
