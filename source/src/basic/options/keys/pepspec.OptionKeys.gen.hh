// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/pepspec.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_pepspec_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_pepspec_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace pepspec { extern BooleanOptionKey const pepspec; }
namespace pepspec { extern StringOptionKey const soft_wts; }
namespace pepspec { extern StringOptionKey const cen_wts; }
namespace pepspec { extern BooleanOptionKey const binding_score; }
namespace pepspec { extern BooleanOptionKey const no_cen; }
namespace pepspec { extern BooleanOptionKey const no_cen_rottrials; }
namespace pepspec { extern BooleanOptionKey const run_sequential; }
namespace pepspec { extern IntegerOptionKey const pep_anchor; }
namespace pepspec { extern StringOptionKey const pep_chain; }
namespace pepspec { extern IntegerOptionKey const n_peptides; }
namespace pepspec { extern IntegerOptionKey const n_build_loop; }
namespace pepspec { extern IntegerOptionKey const n_cgrelax_loop; }
namespace pepspec { extern IntegerOptionKey const n_dock_loop; }
namespace pepspec { extern RealOptionKey const interface_cutoff; }
namespace pepspec { extern BooleanOptionKey const use_input_bb; }
namespace pepspec { extern BooleanOptionKey const remove_input_bb; }
namespace pepspec { extern StringOptionKey const homol_csts; }
namespace pepspec { extern RealOptionKey const p_homol_csts; }
namespace pepspec { extern StringOptionKey const frag_file; }
namespace pepspec { extern BooleanOptionKey const gen_pep_bb_sequential; }
namespace pepspec { extern StringOptionKey const input_seq; }
namespace pepspec { extern StringOptionKey const ss_type; }
namespace pepspec { extern BooleanOptionKey const upweight_interface; }
namespace pepspec { extern BooleanOptionKey const calc_sasa; }
namespace pepspec { extern BooleanOptionKey const diversify_pep_seqs; }
namespace pepspec { extern IntegerOptionKey const diversify_lvl; }
namespace pepspec { extern BooleanOptionKey const dump_cg_bb; }
namespace pepspec { extern BooleanOptionKey const save_low_pdbs; }
namespace pepspec { extern BooleanOptionKey const save_all_pdbs; }
namespace pepspec { extern BooleanOptionKey const no_design; }
namespace pepspec { extern StringOptionKey const pdb_list; }
namespace pepspec { extern StringOptionKey const ref_pdb_list; }
namespace pepspec { extern BooleanOptionKey const add_buffer_res; }
namespace pepspec { extern StringOptionKey const cg_res_type; }
namespace pepspec { extern IntegerOptionKey const native_pep_anchor; }
namespace pepspec { extern StringOptionKey const native_pep_chain; }
namespace pepspec { extern BooleanOptionKey const native_align; }
namespace pepspec { extern BooleanOptionKey const rmsd_analysis; }
namespace pepspec { extern BooleanOptionKey const phipsi_analysis; }
namespace pepspec { extern StringOptionKey const anchor_type; }
namespace pepspec { extern BooleanOptionKey const no_prepack_prot; }
namespace pepspec { extern BooleanOptionKey const prep_use_ref_rotamers; }
namespace pepspec { extern IntegerOptionKey const n_prepend; }
namespace pepspec { extern IntegerOptionKey const n_append; }
namespace pepspec { extern RealOptionKey const clash_cutoff; }
namespace pepspec { extern RealOptionKey const n_anchor_dock_std_devs; }
namespace pepspec { extern RealOptionKey const prep_trans_std_dev; }
namespace pepspec { extern RealOptionKey const prep_rot_std_dev; }
namespace pepspec { extern BooleanOptionKey const seq_align; }
namespace pepspec { extern StringOptionKey const prep_align_prot_to; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
