// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/templates.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_templates_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_templates_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace templates { extern BooleanOptionKey const templates; }
namespace templates { extern FileOptionKey const config; }
namespace templates { extern BooleanOptionKey const fix_aligned_residues; }
namespace templates { extern FileOptionKey const fix_frag_file; }
namespace templates { extern IntegerOptionKey const fix_margin; }
namespace templates { extern IntegerOptionKey const min_nr_large_frags; }
namespace templates { extern IntegerOptionKey const min_nr_small_frags; }
namespace templates { extern BooleanOptionKey const no_pick_fragments; }
namespace templates { extern IntegerOptionKey const nr_large_copies; }
namespace templates { extern IntegerOptionKey const nr_small_copies; }
namespace templates { extern BooleanOptionKey const pairings; }
namespace templates { extern BooleanOptionKey const pick_multiple_sizes; }
namespace templates { extern BooleanOptionKey const strand_constraint; }
namespace templates { extern BooleanOptionKey const vary_frag_size; }
namespace templates { extern BooleanOptionKey const no_culling; }
namespace templates { extern FileOptionKey const helix_pairings; }
namespace templates { extern FileOptionKey const prefix; }
namespace templates { extern IntegerOptionKey const change_movemap; }
namespace templates { extern BooleanOptionKey const force_native_topology; }
namespace templates { extern RealOptionKey const topology_rank_cutoff; }
namespace templates { extern IntegerOptionKey const min_frag_size; }
namespace templates { extern IntegerOptionKey const max_shrink; }
namespace templates { extern IntegerOptionKey const shrink_step; }
namespace templates { extern IntegerOptionKey const shrink_pos_step; }
namespace templates { extern IntegerOptionKey const min_padding; }
namespace templates { extern IntegerOptionKey const min_align_pos; }
namespace templates { extern IntegerOptionKey const max_align_pos; }
namespace templates { namespace cst { extern BooleanOptionKey const cst; } }
namespace templates { namespace cst { extern IntegerOptionKey const topN; } }
namespace templates { namespace cst { extern RealOptionKey const wTopol; } }
namespace templates { namespace cst { extern RealOptionKey const wExtern; } }
namespace templates { namespace fragsteal { extern BooleanOptionKey const fragsteal; } }
namespace templates { namespace fragsteal { extern IntegerOptionKey const topN; } }
namespace templates { namespace fragsteal { extern RealOptionKey const wTopol; } }
namespace templates { namespace fragsteal { extern RealOptionKey const wExtern; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
