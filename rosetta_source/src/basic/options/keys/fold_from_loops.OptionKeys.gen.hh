// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/fold_from_loops.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_fold_from_loops_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_fold_from_loops_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace fold_from_loops { extern BooleanOptionKey const fold_from_loops; }
namespace fold_from_loops { extern BooleanOptionKey const native_ca_cst; }
namespace fold_from_loops { extern FileOptionKey const swap_loops; }
namespace fold_from_loops { extern StringOptionKey const checkpoint; }
namespace fold_from_loops { extern RealOptionKey const ca_csts_dev; }
namespace fold_from_loops { extern IntegerOptionKey const add_relax_cycles; }
namespace fold_from_loops { extern IntegerOptionKey const loop_mov_nterm; }
namespace fold_from_loops { extern IntegerOptionKey const loop_mov_cterm; }
namespace fold_from_loops { extern RealOptionKey const ca_rmsd_cutoff; }
namespace fold_from_loops { extern IntegerVectorOptionKey const res_design_bs; }
namespace fold_from_loops { extern FileOptionKey const clear_csts; }
namespace fold_from_loops { extern BooleanOptionKey const output_centroid; }
namespace fold_from_loops { extern BooleanOptionKey const add_cst_loop; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
