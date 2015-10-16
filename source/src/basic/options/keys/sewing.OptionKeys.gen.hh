// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/sewing.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_sewing_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_sewing_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace sewing { extern BooleanOptionKey const sewing; }
namespace sewing { extern FileOptionKey const model_file_name; }
namespace sewing { extern FileOptionKey const score_file_name; }
namespace sewing { extern FileOptionKey const new_model_file_name; }
namespace sewing { extern StringOptionKey const remove_any_dssp; }
namespace sewing { extern StringOptionKey const remove_all_dssp; }
namespace sewing { extern IntegerOptionKey const min_helix_length; }
namespace sewing { extern IntegerOptionKey const max_helix_length; }
namespace sewing { extern IntegerOptionKey const min_loop_length; }
namespace sewing { extern IntegerOptionKey const max_loop_length; }
namespace sewing { extern IntegerOptionKey const min_strand_length; }
namespace sewing { extern IntegerOptionKey const max_strand_length; }
namespace sewing { extern StringOptionKey const mode; }
namespace sewing { extern IntegerOptionKey const num_models_to_dump; }
namespace sewing { extern IntegerVectorOptionKey const models_to_dump; }
namespace sewing { extern IntegerOptionKey const min_hash_score; }
namespace sewing { extern IntegerOptionKey const max_clash_score; }
namespace sewing { extern IntegerOptionKey const num_segments_to_match; }
namespace sewing { extern IntegerVectorOptionKey const match_segments; }
namespace sewing { extern IntegerOptionKey const max_models; }
namespace sewing { extern IntegerOptionKey const starting_model; }
namespace sewing { extern IntegerOptionKey const num_procs; }
namespace sewing { extern IntegerOptionKey const rank; }
namespace sewing { extern StringOptionKey const assembly_type; }
namespace sewing { extern IntegerOptionKey const num_edges_to_follow; }
namespace sewing { extern RealOptionKey const base_native_bonus; }
namespace sewing { extern IntegerOptionKey const neighbor_cutoff; }
namespace sewing { extern BooleanOptionKey const dump_pdbs; }
namespace sewing { extern BooleanOptionKey const skip_refinement; }
namespace sewing { extern BooleanOptionKey const skip_filters; }
namespace sewing { extern RealOptionKey const min_motif_score; }
namespace sewing { extern IntegerOptionKey const num_repeats; }
namespace sewing { extern BooleanOptionKey const repeat; }
namespace sewing { extern IntegerVectorOptionKey const pose_segment_starts; }
namespace sewing { extern IntegerVectorOptionKey const pose_segment_ends; }
namespace sewing { extern BooleanOptionKey const keep_source_segments; }
namespace sewing { extern FileOptionKey const partner_pdb; }
namespace sewing { extern IntegerVectorOptionKey const keep_model_residues; }
namespace sewing { extern IntegerOptionKey const min_lh_fragments; }
namespace sewing { extern BooleanOptionKey const skip_loop_generation; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
