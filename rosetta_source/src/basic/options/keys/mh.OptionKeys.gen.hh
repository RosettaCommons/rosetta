// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/mh.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_mh_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_mh_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace mh { extern BooleanOptionKey const mh; }
namespace mh { extern StringOptionKey const motif_out_file; }
namespace mh { extern FileVectorOptionKey const harvest_motifs; }
namespace mh { extern FileVectorOptionKey const print_motifs; }
namespace mh { extern FileVectorOptionKey const dump_motif_pdbs; }
namespace mh { extern FileVectorOptionKey const merge_motifs; }
namespace mh { extern BooleanOptionKey const merge_motifs_one_per_bin; }
namespace mh { extern BooleanOptionKey const generate_reverse_motifs; }
namespace mh { extern FileVectorOptionKey const dump_input_pdb; }
namespace mh { extern FileVectorOptionKey const score_pdbs; }
namespace mh { extern FileVectorOptionKey const xform_score_data; }
namespace mh { extern FileVectorOptionKey const xform_score_data_ee; }
namespace mh { extern FileVectorOptionKey const xform_score_data_eh; }
namespace mh { extern FileVectorOptionKey const xform_score_data_he; }
namespace mh { extern FileVectorOptionKey const xform_score_data_hh; }
namespace mh { extern FileVectorOptionKey const xform_score_data_sspair; }
namespace mh { extern FileVectorOptionKey const sequence_recovery; }
namespace mh { extern FileVectorOptionKey const explicit_motif_score; }
namespace mh { extern FileVectorOptionKey const input_motifs; }
namespace mh { extern FileVectorOptionKey const harvest_scores; }
namespace mh { extern FileOptionKey const print_scores; }
namespace mh { extern FileVectorOptionKey const dump_matching_motifs; }
namespace mh { extern RealOptionKey const dump_matching_motifs_cutoff; }
namespace mh { extern BooleanOptionKey const score_across_chains_only; }
namespace mh { extern BooleanOptionKey const normalize_score_ncontact; }
namespace mh { extern IntegerOptionKey const dump_motif_pdbs_min_counts; }
namespace mh { extern RealOptionKey const hash_cart_size; }
namespace mh { extern RealOptionKey const hash_cart_resl; }
namespace mh { extern RealOptionKey const hash_angle_resl; }
namespace mh { extern IntegerOptionKey const harvest_motifs_min_hh_ends; }
namespace mh { extern IntegerOptionKey const harvest_scores_min_count; }
namespace mh { extern BooleanOptionKey const ignore_io_errors; }
namespace mh { extern RealOptionKey const motif_match_radius; }
namespace mh { extern RealVectorOptionKey const merge_similar_motifs; }
namespace mh { namespace score { extern BooleanOptionKey const score; } }
namespace mh { namespace score { extern BooleanOptionKey const noloops; } }
namespace mh { namespace score { extern BooleanOptionKey const spread_ss_element; } }
namespace mh { namespace score { extern RealOptionKey const min_cover_fraction; } }
namespace mh { namespace score { extern RealOptionKey const strand_pair_weight; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter_harvest; } }
namespace mh { namespace filter { extern BooleanOptionKey const filter_io; } }
namespace mh { namespace filter { extern StringOptionKey const restype; } }
namespace mh { namespace filter { extern StringOptionKey const restype_one; } }
namespace mh { namespace filter { extern StringOptionKey const not_restype; } }
namespace mh { namespace filter { extern StringOptionKey const not_restype_one; } }
namespace mh { namespace filter { extern IntegerOptionKey const seqsep; } }
namespace mh { namespace filter { extern BooleanOptionKey const no_hb_bb; } }
namespace mh { namespace filter { extern RealOptionKey const mindist2; } }
namespace mh { namespace filter { extern RealOptionKey const maxdist2; } }
namespace mh { namespace filter { extern StringOptionKey const ss1; } }
namespace mh { namespace filter { extern StringOptionKey const ss2; } }
namespace mh { namespace filter { extern StringOptionKey const dssp1; } }
namespace mh { namespace filter { extern StringOptionKey const dssp2; } }
namespace mh { namespace filter { extern StringOptionKey const aa1; } }
namespace mh { namespace filter { extern StringOptionKey const aa2; } }
namespace mh { namespace filter { extern RealOptionKey const sasa; } }
namespace mh { namespace filter { extern RealOptionKey const faatr; } }
namespace mh { namespace filter { extern RealOptionKey const hb_sc; } }
namespace mh { namespace filter { extern RealOptionKey const hb_bb_sc; } }
namespace mh { namespace filter { extern RealOptionKey const hb_bb; } }
namespace mh { namespace filter { extern RealOptionKey const occupancy; } }
namespace mh { namespace filter { extern RealOptionKey const coorderr; } }
namespace mh { namespace filter { extern RealOptionKey const faatr_or_hbbb; } }
namespace mh { namespace filter { extern RealOptionKey const faatr_or_hb; } }
namespace mh { namespace filter { extern BooleanOptionKey const noloops; } }
namespace mh { namespace filter { extern BooleanOptionKey const oneloop; } }
namespace mh { namespace filter { extern BooleanOptionKey const nodisulf; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
