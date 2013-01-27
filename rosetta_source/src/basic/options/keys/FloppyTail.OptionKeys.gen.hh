// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/FloppyTail.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_FloppyTail_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_FloppyTail_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace FloppyTail { extern BooleanOptionKey const FloppyTail; }
namespace FloppyTail { extern IntegerOptionKey const flexible_start_resnum; }
namespace FloppyTail { extern IntegerOptionKey const flexible_stop_resnum; }
namespace FloppyTail { extern StringOptionKey const flexible_chain; }
namespace FloppyTail { extern RealOptionKey const shear_on; }
namespace FloppyTail { namespace short_tail { extern BooleanOptionKey const short_tail; } }
namespace FloppyTail { namespace short_tail { extern RealOptionKey const short_tail_fraction; } }
namespace FloppyTail { namespace short_tail { extern RealOptionKey const short_tail_off; } }
namespace FloppyTail { extern BooleanOptionKey const pair_off; }
namespace FloppyTail { extern BooleanOptionKey const publication; }
namespace FloppyTail { extern BooleanOptionKey const C_root; }
namespace FloppyTail { extern BooleanOptionKey const force_linear_fold_tree; }
namespace FloppyTail { extern BooleanOptionKey const debug; }
namespace FloppyTail { extern StringOptionKey const cen_weights; }
namespace FloppyTail { extern BooleanOptionKey const perturb_show; }
namespace FloppyTail { extern IntegerOptionKey const perturb_cycles; }
namespace FloppyTail { extern RealOptionKey const perturb_temp; }
namespace FloppyTail { extern IntegerOptionKey const refine_cycles; }
namespace FloppyTail { extern RealOptionKey const refine_temp; }
namespace FloppyTail { extern IntegerOptionKey const refine_repack_cycles; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
