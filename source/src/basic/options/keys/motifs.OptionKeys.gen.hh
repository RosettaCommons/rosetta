// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/motifs.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_motifs_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_motifs_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace motifs { extern BooleanOptionKey const motifs; }
namespace motifs { extern RealOptionKey const close_enough; }
namespace motifs { extern IntegerOptionKey const max_depth; }
namespace motifs { extern BooleanOptionKey const keep_motif_xtal_location; }
namespace motifs { extern RealOptionKey const pack_score_cutoff; }
namespace motifs { extern RealOptionKey const hb_score_cutoff; }
namespace motifs { extern RealOptionKey const water_score_cutoff; }
namespace motifs { extern RealOptionKey const pack_min_threshold; }
namespace motifs { extern RealOptionKey const pack_max_threshold; }
namespace motifs { extern RealOptionKey const hbond_min_threshold; }
namespace motifs { extern RealOptionKey const hbond_max_threshold; }
namespace motifs { extern RealOptionKey const elec_min_threshold; }
namespace motifs { extern RealOptionKey const elec_max_threshold; }
namespace motifs { extern RealOptionKey const duplicate_dist_cutoff; }
namespace motifs { extern RealOptionKey const duplicate_angle_cutoff; }
namespace motifs { extern StringOptionKey const motif_output_directory; }
namespace motifs { extern BooleanOptionKey const eliminate_weak_motifs; }
namespace motifs { extern RealOptionKey const duplicate_motif_cutoff; }
namespace motifs { extern BooleanOptionKey const preminimize_motif_pdbs; }
namespace motifs { extern BooleanOptionKey const preminimize_motif_pdbs_sconly; }
namespace motifs { extern BooleanOptionKey const place_adduct_waters; }
namespace motifs { extern FileVectorOptionKey const list_motifs; }
namespace motifs { extern StringOptionKey const motif_filename; }
namespace motifs { extern StringOptionKey const BPData; }
namespace motifs { extern FileVectorOptionKey const list_dnaconformers; }
namespace motifs { extern StringVectorOptionKey const target_dna_defs; }
namespace motifs { extern StringVectorOptionKey const motif_build_defs; }
namespace motifs { extern StringOptionKey const motif_build_position_chain; }
namespace motifs { extern IntegerVectorOptionKey const motif_build_positions; }
namespace motifs { extern RealOptionKey const r1; }
namespace motifs { extern RealOptionKey const r2; }
namespace motifs { extern RealOptionKey const z1; }
namespace motifs { extern RealOptionKey const z2; }
namespace motifs { extern RealOptionKey const dtest; }
namespace motifs { extern IntegerOptionKey const rotlevel; }
namespace motifs { extern IntegerOptionKey const num_repacks; }
namespace motifs { extern BooleanOptionKey const minimize; }
namespace motifs { extern BooleanOptionKey const minimize_dna; }
namespace motifs { extern BooleanOptionKey const run_motifs; }
namespace motifs { extern BooleanOptionKey const expand_motifs; }
namespace motifs { extern BooleanOptionKey const aromatic_motifs; }
namespace motifs { extern BooleanOptionKey const dump_motifs; }
namespace motifs { extern BooleanOptionKey const quick_and_dirty; }
namespace motifs { extern RealOptionKey const special_rotweight; }
namespace motifs { extern StringOptionKey const output_file; }
namespace motifs { extern StringOptionKey const data_file; }
namespace motifs { extern StringOptionKey const target_aa; }
namespace motifs { extern RealOptionKey const constraint_max; }
namespace motifs { extern BooleanOptionKey const flex_sugar; }
namespace motifs { extern BooleanOptionKey const clear_bprots; }
namespace motifs { extern IntegerOptionKey const rots2add; }
namespace motifs { extern BooleanOptionKey const restrict_to_wt; }
namespace motifs { extern BooleanOptionKey const rerun_motifsearch; }
namespace motifs { extern BooleanOptionKey const no_rotamer_bump; }
namespace motifs { extern RealOptionKey const ligand_motif_sphere; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
