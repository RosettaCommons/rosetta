// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/jumps.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_jumps_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_jumps_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace jumps { extern BooleanOptionKey const jumps; }
namespace jumps { extern BooleanOptionKey const evaluate; }
namespace jumps { extern FileOptionKey const extra_frags_for_ss; }
namespace jumps { extern BooleanOptionKey const fix_chainbreak; }
namespace jumps { extern FileOptionKey const fix_jumps; }
namespace jumps { extern FileOptionKey const jump_lib; }
namespace jumps { extern FileOptionKey const loop_definition_from_file; }
namespace jumps { extern BooleanOptionKey const no_chainbreak_in_relax; }
namespace jumps { extern FileOptionKey const pairing_file; }
namespace jumps { extern IntegerVectorOptionKey const random_sheets; }
namespace jumps { extern FileOptionKey const residue_pair_jump_file; }
namespace jumps { extern IntegerVectorOptionKey const sheets; }
namespace jumps { extern FileOptionKey const topology_file; }
namespace jumps { extern BooleanOptionKey const bb_moves; }
namespace jumps { extern BooleanOptionKey const no_wobble; }
namespace jumps { extern BooleanOptionKey const no_shear; }
namespace jumps { extern BooleanOptionKey const no_sample_ss_jumps; }
namespace jumps { extern IntegerOptionKey const invrate_jump_move; }
namespace jumps { extern RealOptionKey const chainbreak_weight_stage1; }
namespace jumps { extern RealOptionKey const chainbreak_weight_stage2; }
namespace jumps { extern RealOptionKey const chainbreak_weight_stage3; }
namespace jumps { extern RealOptionKey const chainbreak_weight_stage4; }
namespace jumps { extern BooleanOptionKey const ramp_chainbreaks; }
namespace jumps { extern RealOptionKey const increase_chainbreak; }
namespace jumps { extern BooleanOptionKey const overlap_chainbreak; }
namespace jumps { extern RealOptionKey const sep_switch_accelerate; }
namespace jumps { extern BooleanOptionKey const dump_frags; }
namespace jumps { extern IntegerOptionKey const njumps; }
namespace jumps { extern IntegerOptionKey const max_strand_gap_allowed; }
namespace jumps { extern RealOptionKey const contact_score; }
namespace jumps { extern BooleanOptionKey const filter_templates; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
