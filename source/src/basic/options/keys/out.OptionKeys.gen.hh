// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/out.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_out_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_out_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace out { extern BooleanOptionKey const out; }
namespace out { extern BooleanOptionKey const overwrite; }
namespace out { extern IntegerOptionKey const nstruct; }
namespace out { extern IntegerOptionKey const shuffle_nstruct; }
namespace out { extern StringOptionKey const prefix; }
namespace out { extern StringOptionKey const suffix; }
namespace out { extern StringOptionKey const force_output_name; }
namespace out { extern BooleanOptionKey const no_nstruct_label; }
namespace out { extern BooleanOptionKey const pdb_gz; }
namespace out { extern BooleanOptionKey const pdb; }
namespace out { extern BooleanOptionKey const silent_gz; }
namespace out { extern BooleanOptionKey const use_database; }
namespace out { extern IntegerOptionKey const database_protocol_id; }
namespace out { extern StringVectorOptionKey const database_filter; }
namespace out { extern IntegerVectorOptionKey const resume_batch; }
namespace out { extern BooleanOptionKey const nooutput; }
namespace out { extern BooleanOptionKey const output; }
namespace out { extern RealOptionKey const scorecut; }
namespace out { extern BooleanOptionKey const show_accessed_options; }
namespace out { extern FileOptionKey const sf; }
namespace out { extern StringVectorOptionKey const mute; }
namespace out { extern StringVectorOptionKey const unmute; }
namespace out { extern IntegerOptionKey const level; }
namespace out { extern StringVectorOptionKey const levels; }
namespace out { extern IntegerOptionKey const std_IO_exit_error_code; }
namespace out { extern BooleanOptionKey const chname; }
namespace out { extern BooleanOptionKey const chtimestamp; }
namespace out { extern BooleanOptionKey const dry_run; }
namespace out { extern StringOptionKey const mpi_tracer_to_file; }
namespace out { extern StringOptionKey const user_tag; }
namespace out { extern StringOptionKey const output_tag; }
namespace out { namespace file { extern BooleanOptionKey const file; } }
namespace out { namespace file { extern StringOptionKey const o; } }
namespace out { namespace file { extern StringOptionKey const silent; } }
namespace out { namespace file { extern StringOptionKey const score_only; } }
namespace out { namespace file { extern StringOptionKey const atom_tree_diff; } }
namespace out { namespace file { extern IntegerOptionKey const atom_tree_diff_bb; } }
namespace out { namespace file { extern IntegerOptionKey const atom_tree_diff_sc; } }
namespace out { namespace file { extern IntegerOptionKey const atom_tree_diff_bl; } }
namespace out { namespace file { extern StringOptionKey const alignment; } }
namespace out { namespace file { extern StringOptionKey const scorefile; } }
namespace out { namespace file { extern StringOptionKey const silent_struct_type; } }
namespace out { namespace file { extern BooleanOptionKey const silent_preserve_H; } }
namespace out { namespace file { extern BooleanOptionKey const silent_print_all_score_headers; } }
namespace out { namespace file { extern BooleanOptionKey const silent_decoytime; } }
namespace out { namespace file { extern IntegerOptionKey const silent_comment_bound; } }
namespace out { namespace file { extern BooleanOptionKey const raw; } }
namespace out { namespace file { extern BooleanOptionKey const fullatom; } }
namespace out { namespace file { extern BooleanOptionKey const suppress_zero_occ_pdb_output; } }
namespace out { namespace file { extern BooleanOptionKey const output_virtual; } }
namespace out { namespace file { extern BooleanOptionKey const no_output_cen; } }
namespace out { namespace file { extern BooleanOptionKey const output_orbitals; } }
namespace out { namespace file { extern BooleanOptionKey const weight_silent_scores; } }
namespace out { namespace file { extern FileOptionKey const design_contrast; } }
namespace out { namespace file { extern BooleanOptionKey const dont_rewrite_dunbrack_database; } }
namespace out { namespace file { extern BooleanOptionKey const renumber_pdb; } }
namespace out { namespace file { extern BooleanOptionKey const pdb_parents; } }
namespace out { namespace file { extern BooleanOptionKey const per_chain_renumbering; } }
namespace out { namespace file { extern StringOptionKey const residue_type_set; } }
namespace out { namespace file { extern StringOptionKey const frag_prefix; } }
namespace out { namespace file { extern BooleanOptionKey const output_torsions; } }
namespace out { namespace file { extern BooleanOptionKey const pdb_comments; } }
namespace out { namespace file { extern BooleanOptionKey const force_nonideal_structure; } }
namespace out { namespace path { extern PathOptionKey const all; } }
namespace out { namespace path { extern PathOptionKey const path; } }
namespace out { namespace path { extern PathOptionKey const pdb; } }
namespace out { namespace path { extern PathOptionKey const score; } }
namespace out { namespace path { extern PathOptionKey const movie; } }
namespace out { namespace path { extern PathOptionKey const scratch; } }
namespace out { namespace path { extern BooleanOptionKey const mpi_rank_dir; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
