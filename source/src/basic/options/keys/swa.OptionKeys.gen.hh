// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/swa.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_swa_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_swa_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace swa { extern BooleanOptionKey const swa; }
namespace swa { extern StringVectorOptionKey const s1; }
namespace swa { extern StringVectorOptionKey const s2; }
namespace swa { extern StringVectorOptionKey const silent1; }
namespace swa { extern StringVectorOptionKey const silent2; }
namespace swa { extern StringVectorOptionKey const tags1; }
namespace swa { extern StringVectorOptionKey const tags2; }
namespace swa { extern IntegerVectorOptionKey const slice_res1; }
namespace swa { extern IntegerVectorOptionKey const slice_res2; }
namespace swa { extern IntegerVectorOptionKey const input_res1; }
namespace swa { extern IntegerVectorOptionKey const input_res2; }
namespace swa { extern BooleanOptionKey const backbone_only1; }
namespace swa { extern BooleanOptionKey const backbone_only2; }
namespace swa { namespace rna { extern BooleanOptionKey const rna; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const sample_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const fixed_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const minimize_res; } }
namespace swa { namespace rna { extern StringVectorOptionKey const alignment_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const native_alignment_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const rmsd_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const missing_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const missing_res2; } }
namespace swa { namespace rna { extern IntegerOptionKey const job_queue_ID; } }
namespace swa { namespace rna { extern BooleanOptionKey const minimize_and_score_sugar; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const global_sample_res_list; } }
namespace swa { namespace rna { extern FileOptionKey const filter_output_filename; } }
namespace swa { namespace rna { extern BooleanOptionKey const combine_long_loop_mode; } }
namespace swa { namespace rna { extern BooleanOptionKey const combine_helical_silent_file; } }
namespace swa { namespace rna { extern BooleanOptionKey const output_extra_RMSDs; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const force_syn_chi_res_list; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const force_north_sugar_list; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const force_south_sugar_list; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const protonated_H1_adenosine_list; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const native_virtual_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const simple_append_map; } }
namespace swa { namespace rna { extern BooleanOptionKey const allow_fixed_res_at_moving_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const force_user_defined_jumps; } }
namespace swa { namespace rna { extern BooleanOptionKey const test_encapsulation; } }
namespace swa { namespace rna { extern StringVectorOptionKey const jump_point_pairs; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const terminal_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const add_virt_root; } }
namespace swa { namespace rna { extern BooleanOptionKey const floating_base; } }
namespace swa { namespace rna { extern IntegerOptionKey const floating_base_anchor_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const allow_chain_boundary_jump_partner_right_at_fixed_BP; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const virtual_res; } }
namespace swa { namespace rna { extern IntegerVectorOptionKey const bulge_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const rebuild_bulge_mode; } }
namespace swa { namespace rna { extern BooleanOptionKey const choose_random; } }
namespace swa { namespace rna { extern IntegerOptionKey const num_random_samples; } }
namespace swa { namespace rna { extern BooleanOptionKey const filter_user_alignment_res; } }
namespace swa { namespace rna { extern BooleanOptionKey const output_pdb; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
