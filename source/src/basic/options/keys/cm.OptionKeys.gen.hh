// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/cm.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_cm_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_cm_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace cm { extern BooleanOptionKey const cm; }
namespace cm { namespace sanitize { extern BooleanOptionKey const sanitize; } }
namespace cm { namespace sanitize { extern RealOptionKey const bound_delta; } }
namespace cm { namespace sanitize { extern RealOptionKey const bound_sd; } }
namespace cm { namespace sanitize { extern IntegerOptionKey const num_fragments; } }
namespace cm { namespace sanitize { extern RealOptionKey const cst_weight_pair; } }
namespace cm { namespace sanitize { extern RealOptionKey const cst_weight_coord; } }
namespace cm { extern BooleanOptionKey const start_models_only; }
namespace cm { extern StringOptionKey const aln_format; }
namespace cm { extern BooleanOptionKey const recover_side_chains; }
namespace cm { extern FileVectorOptionKey const steal_extra_residues; }
namespace cm { extern StringOptionKey const loop_mover; }
namespace cm { extern IntegerOptionKey const loop_close_level; }
namespace cm { extern IntegerOptionKey const min_loop_size; }
namespace cm { extern IntegerOptionKey const max_loop_rebuild; }
namespace cm { extern RealOptionKey const loop_rebuild_filter; }
namespace cm { extern RealOptionKey const aln_length_filter_quantile; }
namespace cm { extern IntegerOptionKey const aln_length_filter; }
namespace cm { extern StringVectorOptionKey const template_ids; }
namespace cm { extern FileOptionKey const ligand_pdb; }
namespace cm { extern StringVectorOptionKey const seq_score; }
namespace cm { extern StringOptionKey const aligner; }
namespace cm { extern RealOptionKey const min_gap_open; }
namespace cm { extern RealOptionKey const max_gap_open; }
namespace cm { extern RealOptionKey const min_gap_extend; }
namespace cm { extern RealOptionKey const max_gap_extend; }
namespace cm { extern IntegerOptionKey const nn; }
namespace cm { extern RealOptionKey const fr_temperature; }
namespace cm { extern FileVectorOptionKey const ev_map; }
namespace cm { extern FileVectorOptionKey const hh_map; }
namespace cm { namespace hybridize { extern BooleanOptionKey const hybridize; } }
namespace cm { namespace hybridize { extern FileVectorOptionKey const templates; } }
namespace cm { namespace hybridize { extern FileOptionKey const template_list; } }
namespace cm { namespace hybridize { extern IntegerVectorOptionKey const starting_template; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const realign_domains; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const add_non_init_chunks; } }
namespace cm { namespace hybridize { extern StringOptionKey const ss; } }
namespace cm { namespace hybridize { extern RealOptionKey const stage1_increase_cycles; } }
namespace cm { namespace hybridize { extern RealOptionKey const stage2_increase_cycles; } }
namespace cm { namespace hybridize { extern RealOptionKey const stage2min_increase_cycles; } }
namespace cm { namespace hybridize { extern RealOptionKey const stage1_probability; } }
namespace cm { namespace hybridize { extern StringOptionKey const stage1_weights; } }
namespace cm { namespace hybridize { extern StringOptionKey const stage1_patch; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const skip_stage2; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const no_global_frame; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const linmin_only; } }
namespace cm { namespace hybridize { extern StringOptionKey const stage2_weights; } }
namespace cm { namespace hybridize { extern StringOptionKey const stage2_patch; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const relax; } }
namespace cm { namespace hybridize { extern RealOptionKey const frag_weight_aligned; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const move_anchor; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const max_registry_shift; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const alignment_from_template_seqpos; } }
namespace cm { namespace hybridize { extern IntegerVectorOptionKey const alignment_from_chunk_mapping; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const virtual_loops; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const revert_real_loops; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const realign_domains_stage2; } }
namespace cm { namespace hybridize { extern RealOptionKey const frag_1mer_insertion_weight; } }
namespace cm { namespace hybridize { extern RealOptionKey const small_frag_insertion_weight; } }
namespace cm { namespace hybridize { extern RealOptionKey const big_frag_insertion_weight; } }
namespace cm { namespace hybridize { extern BooleanOptionKey const auto_frag_insertion_weight; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const stage1_1_cycles; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const stage1_2_cycles; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const stage1_3_cycles; } }
namespace cm { namespace hybridize { extern IntegerOptionKey const stage1_4_cycles; } }
namespace cm { namespace hybridize { extern RealOptionKey const stage2_temperature; } }
namespace cm { namespace hybridize { extern StringOptionKey const stage1_4_cenrot_score; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
