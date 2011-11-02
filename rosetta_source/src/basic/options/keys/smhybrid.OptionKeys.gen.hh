// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/smhybrid.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_smhybrid_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_smhybrid_OptionKeys_gen_HH

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

namespace smhybrid { extern BooleanOptionKey const smhybrid; }
namespace smhybrid { extern BooleanOptionKey const add_cavities; }
namespace smhybrid { extern BooleanOptionKey const abinitio_design; }
namespace smhybrid { extern BooleanOptionKey const fa_refine; }
namespace smhybrid { extern BooleanOptionKey const virtual_nterm; }
namespace smhybrid { extern BooleanOptionKey const debug; }
namespace smhybrid { extern BooleanOptionKey const refine; }
namespace smhybrid { extern BooleanOptionKey const filter; }
namespace smhybrid { extern BooleanOptionKey const floating_scs_rep; }
namespace smhybrid { extern BooleanOptionKey const flxbb; }
namespace smhybrid { extern BooleanOptionKey const centroid_all_val; }
namespace smhybrid { extern BooleanOptionKey const subsubs_attract; }
namespace smhybrid { extern BooleanOptionKey const linker_cst; }
namespace smhybrid { extern BooleanOptionKey const pseudosym; }
namespace smhybrid { extern BooleanOptionKey const design_linker; }
namespace smhybrid { extern BooleanOptionKey const design; }
namespace smhybrid { extern BooleanOptionKey const restrict_design_to_interface; }
namespace smhybrid { extern BooleanOptionKey const restrict_design_to_subsub_interface; }
namespace smhybrid { extern BooleanOptionKey const design_hydrophobic; }
namespace smhybrid { extern BooleanOptionKey const add_metal_at_0; }
namespace smhybrid { extern IntegerOptionKey const nres_mono; }
namespace smhybrid { extern IntegerOptionKey const abinitio_cycles; }
namespace smhybrid { extern IntegerOptionKey const primary_subsubunit; }
namespace smhybrid { extern IntegerOptionKey const minbb; }
namespace smhybrid { extern IntegerOptionKey const switch_concert_sub; }
namespace smhybrid { extern RealOptionKey const temperature; }
namespace smhybrid { extern BooleanOptionKey const inter_subsub_cst; }
namespace smhybrid { extern RealOptionKey const rb_mag; }
namespace smhybrid { extern StringOptionKey const ss; }
namespace smhybrid { extern FileOptionKey const symm_def_template; }
namespace smhybrid { extern FileOptionKey const symm_def_template_reduced; }
namespace smhybrid { extern IntegerVectorOptionKey const attach_as_sc; }
namespace smhybrid { extern IntegerVectorOptionKey const attach_as_sc_sub; }
namespace smhybrid { extern IntegerVectorOptionKey const inversion_subs; }
namespace smhybrid { extern BooleanVectorOptionKey const chainbreaks; }
namespace smhybrid { extern StringVectorOptionKey const design_res_files; }
namespace smhybrid { extern StringVectorOptionKey const fixed_res_files; }
namespace smhybrid { extern StringVectorOptionKey const frag_res_files; }
namespace smhybrid { extern StringVectorOptionKey const scattach_res_files; }
namespace smhybrid { extern StringVectorOptionKey const rep_edge_files; }
namespace smhybrid { extern StringVectorOptionKey const virtual_res_files; }
namespace smhybrid { extern StringVectorOptionKey const jumpcut_files; }
namespace smhybrid { extern StringVectorOptionKey const cst_sub_files; }
namespace smhybrid { extern StringVectorOptionKey const symm_file_tag; }
namespace smhybrid { extern StringVectorOptionKey const attach_atom; }
namespace smhybrid { extern StringVectorOptionKey const add_res_before; }
namespace smhybrid { extern StringVectorOptionKey const add_res_after; }
namespace smhybrid { extern StringVectorOptionKey const add_ss_before; }
namespace smhybrid { extern StringVectorOptionKey const add_ss_after; }
namespace smhybrid { extern StringVectorOptionKey const add_atom_at_cen; }
namespace smhybrid { extern StringVectorOptionKey const attach_rsd; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
