// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/flexPepDocking.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_flexPepDocking_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_flexPepDocking_OptionKeys_gen_HH

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

namespace flexPepDocking { extern BooleanOptionKey const flexPepDocking; }
namespace flexPepDocking { extern StringOptionKey const params_file; }
namespace flexPepDocking { extern IntegerOptionKey const peptide_anchor; }
namespace flexPepDocking { extern StringOptionKey const receptor_chain; }
namespace flexPepDocking { extern StringOptionKey const peptide_chain; }
namespace flexPepDocking { extern BooleanOptionKey const pep_fold_only; }
namespace flexPepDocking { extern BooleanOptionKey const lowres_abinitio; }
namespace flexPepDocking { extern BooleanOptionKey const lowres_preoptimize; }
namespace flexPepDocking { extern BooleanOptionKey const flexPepDockingMinimizeOnly; }
namespace flexPepDocking { extern BooleanOptionKey const extend_peptide; }
namespace flexPepDocking { extern BooleanOptionKey const pep_refine; }
namespace flexPepDocking { extern BooleanOptionKey const rbMCM; }
namespace flexPepDocking { extern BooleanOptionKey const torsionsMCM; }
namespace flexPepDocking { extern BooleanOptionKey const peptide_loop_model; }
namespace flexPepDocking { extern BooleanOptionKey const backrub_peptide; }
namespace flexPepDocking { extern BooleanOptionKey const boost_fa_atr; }
namespace flexPepDocking { extern BooleanOptionKey const ramp_fa_rep; }
namespace flexPepDocking { extern BooleanOptionKey const ramp_rama; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_score_only; }
namespace flexPepDocking { extern FileOptionKey const ref_startstruct; }
namespace flexPepDocking { extern BooleanOptionKey const use_cen_score; }
namespace flexPepDocking { extern BooleanOptionKey const design_peptide; }
namespace flexPepDocking { extern IntegerOptionKey const rep_ramp_cycles; }
namespace flexPepDocking { extern IntegerOptionKey const mcm_cycles; }
namespace flexPepDocking { extern RealOptionKey const random_phi_psi_preturbation; }
namespace flexPepDocking { extern RealOptionKey const smove_angle_range; }
namespace flexPepDocking { extern BooleanOptionKey const min_receptor_bb; }
namespace flexPepDocking { extern RealOptionKey const random_trans_start; }
namespace flexPepDocking { extern RealOptionKey const random_rot_start; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_prepack; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_noprepack1; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_noprepack2; }
namespace flexPepDocking { extern RealOptionKey const score_filter; }
namespace flexPepDocking { extern IntegerOptionKey const hb_filter; }
namespace flexPepDocking { extern IntegerOptionKey const hotspot_filter; }
namespace flexPepDocking { extern StringOptionKey const frag5; }
namespace flexPepDocking { extern RealOptionKey const frag9_weight; }
namespace flexPepDocking { extern RealOptionKey const frag5_weight; }
namespace flexPepDocking { extern RealOptionKey const frag3_weight; }
namespace flexPepDocking { extern BooleanOptionKey const pSer2Asp_centroid; }
namespace flexPepDocking { extern BooleanOptionKey const pSer2Glu_centroid; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
