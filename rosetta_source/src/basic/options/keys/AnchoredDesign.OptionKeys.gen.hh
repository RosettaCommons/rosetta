// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/AnchoredDesign.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_AnchoredDesign_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_AnchoredDesign_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace AnchoredDesign { extern BooleanOptionKey const AnchoredDesign; }
namespace AnchoredDesign { extern FileOptionKey const anchor; }
namespace AnchoredDesign { extern BooleanOptionKey const allow_anchor_repack; }
namespace AnchoredDesign { extern BooleanOptionKey const vary_cutpoints; }
namespace AnchoredDesign { extern BooleanOptionKey const no_frags; }
namespace AnchoredDesign { extern BooleanOptionKey const debug; }
namespace AnchoredDesign { extern BooleanOptionKey const show_extended; }
namespace AnchoredDesign { extern BooleanOptionKey const refine_only; }
namespace AnchoredDesign { extern BooleanOptionKey const perturb_show; }
namespace AnchoredDesign { extern IntegerOptionKey const perturb_cycles; }
namespace AnchoredDesign { extern RealOptionKey const perturb_temp; }
namespace AnchoredDesign { extern BooleanOptionKey const perturb_CCD_off; }
namespace AnchoredDesign { extern BooleanOptionKey const perturb_KIC_off; }
namespace AnchoredDesign { extern BooleanOptionKey const refine_CCD_off; }
namespace AnchoredDesign { extern BooleanOptionKey const refine_KIC_off; }
namespace AnchoredDesign { extern IntegerOptionKey const refine_cycles; }
namespace AnchoredDesign { extern RealOptionKey const refine_temp; }
namespace AnchoredDesign { extern IntegerOptionKey const refine_repack_cycles; }
namespace AnchoredDesign { extern BooleanOptionKey const rmsd; }
namespace AnchoredDesign { extern BooleanOptionKey const unbound_mode; }
namespace AnchoredDesign { extern RealOptionKey const chainbreak_weight; }
namespace AnchoredDesign { namespace filters { extern BooleanOptionKey const filters; } }
namespace AnchoredDesign { namespace filters { extern RealOptionKey const score; } }
namespace AnchoredDesign { namespace filters { extern RealOptionKey const sasa; } }
namespace AnchoredDesign { namespace filters { extern BooleanOptionKey const omega; } }
namespace AnchoredDesign { namespace akash { extern BooleanOptionKey const akash; } }
namespace AnchoredDesign { namespace akash { extern IntegerOptionKey const dyepos; } }
namespace AnchoredDesign { namespace testing { extern BooleanOptionKey const testing; } }
namespace AnchoredDesign { namespace testing { extern RealOptionKey const VDW_weight; } }
namespace AnchoredDesign { namespace testing { extern BooleanOptionKey const anchor_via_constraints; } }
namespace AnchoredDesign { namespace testing { extern BooleanOptionKey const delete_interface_native_sidechains; } }
namespace AnchoredDesign { namespace testing { extern FileOptionKey const RMSD_only_this; } }
namespace AnchoredDesign { namespace testing { extern BooleanOptionKey const anchor_noise_constraints_mode; } }
namespace AnchoredDesign { namespace testing { extern BooleanOptionKey const super_secret_fixed_interface_mode; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
