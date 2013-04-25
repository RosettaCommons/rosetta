// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/enzdes.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_enzdes_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_enzdes_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace enzdes { extern BooleanOptionKey const enzdes; }
namespace enzdes { extern StringOptionKey const checkpoint; }
namespace enzdes { extern BooleanOptionKey const enz_score; }
namespace enzdes { extern BooleanOptionKey const enz_repack; }
namespace enzdes { extern BooleanOptionKey const cst_opt; }
namespace enzdes { extern BooleanOptionKey const cst_predock; }
namespace enzdes { extern RealOptionKey const trans_magnitude; }
namespace enzdes { extern RealOptionKey const rot_magnitude; }
namespace enzdes { extern RealOptionKey const dock_trials; }
namespace enzdes { extern BooleanOptionKey const cst_min; }
namespace enzdes { extern BooleanOptionKey const cst_design; }
namespace enzdes { extern IntegerOptionKey const design_min_cycles; }
namespace enzdes { extern BooleanOptionKey const make_consensus_mutations; }
namespace enzdes { extern BooleanOptionKey const bb_min; }
namespace enzdes { extern RealOptionKey const bb_min_allowed_dev; }
namespace enzdes { extern RealOptionKey const loop_bb_min_allowed_dev; }
namespace enzdes { extern RealOptionKey const minimize_ligand_torsions; }
namespace enzdes { extern RealOptionKey const minimize_all_ligand_torsions; }
namespace enzdes { extern BooleanOptionKey const chi_min; }
namespace enzdes { extern BooleanOptionKey const min_all_jumps; }
namespace enzdes { extern BooleanOptionKey const cst_dock; }
namespace enzdes { extern BooleanOptionKey const run_ligand_motifs; }
namespace enzdes { extern BooleanOptionKey const enz_debug; }
namespace enzdes { extern FileOptionKey const cstfile; }
namespace enzdes { extern FileOptionKey const enz_loops_file; }
namespace enzdes { extern BooleanOptionKey const flexbb_protocol; }
namespace enzdes { extern BooleanOptionKey const remodel_protocol; }
namespace enzdes { extern BooleanOptionKey const kic_loop_sampling; }
namespace enzdes { extern StringOptionKey const dump_loop_samples; }
namespace enzdes { extern BooleanOptionKey const fix_catalytic_aa; }
namespace enzdes { extern IntegerOptionKey const additional_packing_ligand_rb_confs; }
namespace enzdes { extern IntegerOptionKey const ex_catalytic_rot; }
namespace enzdes { extern IntegerOptionKey const single_loop_ensemble_size; }
namespace enzdes { extern IntegerOptionKey const loop_generator_trials; }
namespace enzdes { extern BooleanOptionKey const no_catres_min_in_loopgen; }
namespace enzdes { extern RealOptionKey const mc_kt_low; }
namespace enzdes { extern RealOptionKey const mc_kt_high; }
namespace enzdes { extern RealOptionKey const min_cacb_deviation; }
namespace enzdes { extern RealOptionKey const max_bb_deviation; }
namespace enzdes { extern RealOptionKey const max_bb_deviation_from_startstruct; }
namespace enzdes { extern IntegerOptionKey const flexbb_outstructs; }
namespace enzdes { extern IntegerOptionKey const remodel_trials; }
namespace enzdes { extern BooleanOptionKey const remodel_secmatch; }
namespace enzdes { extern BooleanOptionKey const dump_inverse_rotamers; }
namespace enzdes { extern RealOptionKey const remodel_aggressiveness; }
namespace enzdes { extern RealOptionKey const favor_native_res; }
namespace enzdes { extern BooleanOptionKey const detect_design_interface; }
namespace enzdes { extern BooleanOptionKey const include_catres_in_interface_detection; }
namespace enzdes { extern BooleanOptionKey const arg_sweep_interface; }
namespace enzdes { extern RealOptionKey const arg_sweep_cutoff; }
namespace enzdes { extern RealOptionKey const cut1; }
namespace enzdes { extern RealOptionKey const cut2; }
namespace enzdes { extern RealOptionKey const cut3; }
namespace enzdes { extern RealOptionKey const cut4; }
namespace enzdes { extern RealOptionKey const lig_packer_weight; }
namespace enzdes { extern BooleanOptionKey const no_unconstrained_repack; }
namespace enzdes { extern RealOptionKey const secmatch_Ecutoff; }
namespace enzdes { extern FileOptionKey const change_lig; }
namespace enzdes { extern StringOptionKey const process_ligrot_separately; }
namespace enzdes { extern BooleanOptionKey const start_from_random_rb_conf; }
namespace enzdes { extern RealOptionKey const bb_bump_cutoff; }
namespace enzdes { extern RealOptionKey const sc_sc_bump_cutoff; }
namespace enzdes { extern BooleanOptionKey const no_packstat_calculation; }
namespace enzdes { extern StringOptionKey const compare_native; }
namespace enzdes { extern BooleanOptionKey const final_repack_without_ligand; }
namespace enzdes { extern BooleanOptionKey const dump_final_repack_without_ligand_pdb; }
namespace enzdes { extern BooleanOptionKey const parser_read_cloud_pdb; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
