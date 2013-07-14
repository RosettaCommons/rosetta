// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/optE.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_optE_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_optE_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace optE { extern BooleanOptionKey const optE; }
namespace optE { extern StringOptionKey const load_from_silent; }
namespace optE { extern StringOptionKey const data_in; }
namespace optE { extern StringOptionKey const data_out; }
namespace optE { extern StringOptionKey const weights; }
namespace optE { extern StringVectorOptionKey const fix; }
namespace optE { extern FileOptionKey const free; }
namespace optE { extern FileOptionKey const fixed; }
namespace optE { extern FileOptionKey const parse_tagfile; }
namespace optE { extern FileOptionKey const constant_logic_taskops_file; }
namespace optE { extern BooleanOptionKey const optE_soft_rep; }
namespace optE { extern BooleanOptionKey const no_hb_env_dependence; }
namespace optE { extern BooleanOptionKey const no_hb_env_dependence_DNA; }
namespace optE { extern BooleanOptionKey const optE_no_protein_fa_elec; }
namespace optE { extern BooleanOptionKey const design_first; }
namespace optE { extern IntegerOptionKey const n_design_cycles; }
namespace optE { extern BooleanOptionKey const recover_nat_rot; }
namespace optE { extern FileOptionKey const component_weights; }
namespace optE { extern BooleanOptionKey const optimize_nat_aa; }
namespace optE { extern BooleanOptionKey const optimize_nat_rot; }
namespace optE { extern FileOptionKey const optimize_ligand_rot; }
namespace optE { extern BooleanOptionKey const optimize_pssm; }
namespace optE { extern FileOptionKey const optimize_dGbinding; }
namespace optE { extern FileOptionKey const optimize_ddG_bind_correlation; }
namespace optE { extern FileOptionKey const optimize_ddGmutation; }
namespace optE { extern BooleanOptionKey const optimize_ddGmutation_straight_mean; }
namespace optE { extern BooleanOptionKey const optimize_ddGmutation_boltzman_average; }
namespace optE { extern RealOptionKey const exclude_badrep_ddGs; }
namespace optE { extern BooleanOptionKey const pretend_no_ddG_repulsion; }
namespace optE { extern FileOptionKey const optimize_decoy_discrimination; }
namespace optE { extern StringOptionKey const normalize_decoy_score_spread; }
namespace optE { extern BooleanOptionKey const ramp_nativeness; }
namespace optE { extern IntegerOptionKey const n_top_natives_to_optimize; }
namespace optE { extern RealOptionKey const approximate_decoy_entropy; }
namespace optE { extern BooleanOptionKey const repack_and_minimize_decoys; }
namespace optE { extern BooleanOptionKey const repack_and_minimize_input_structures; }
namespace optE { extern IntegerOptionKey const output_top_n_new_decoys; }
namespace optE { extern FileOptionKey const optimize_ligand_discrimination; }
namespace optE { extern BooleanOptionKey const no_design; }
namespace optE { extern BooleanOptionKey const sqrt_pssm; }
namespace optE { extern RealOptionKey const min_decoy_rms_to_native; }
namespace optE { extern RealOptionKey const max_rms_from_native; }
namespace optE { extern BooleanOptionKey const optimize_starting_free_weights; }
namespace optE { extern FileOptionKey const wrap_dof_optimization; }
namespace optE { extern RealOptionKey const randomly_perturb_starting_free_weights; }
namespace optE { extern RealOptionKey const inv_kT_natrot; }
namespace optE { extern RealOptionKey const inv_kT_nataa; }
namespace optE { extern RealOptionKey const inv_kT_natstruct; }
namespace optE { extern BooleanOptionKey const mpi_weight_minimization; }
namespace optE { extern BooleanOptionKey const dont_use_reference_energies; }
namespace optE { extern IntegerOptionKey const number_of_swarm_particles; }
namespace optE { extern IntegerOptionKey const number_of_swarm_cycles; }
namespace optE { extern FileOptionKey const constrain_weights; }
namespace optE { extern BooleanOptionKey const fit_reference_energies_to_aa_profile_recovery; }
namespace optE { extern FileOptionKey const starting_refEs; }
namespace optE { extern BooleanOptionKey const repeat_swarm_optimization_until_fitness_improves; }
namespace optE { extern BooleanOptionKey const design_with_minpack; }
namespace optE { extern BooleanOptionKey const limit_bad_scores; }
namespace optE { namespace rescore { extern BooleanOptionKey const rescore; } }
namespace optE { namespace rescore { extern FileOptionKey const weights; } }
namespace optE { namespace rescore { extern IntegerOptionKey const context_round; } }
namespace optE { namespace rescore { extern FileOptionKey const outlog; } }
namespace optE { namespace rescore { extern BooleanOptionKey const measure_sequence_recovery; } }
namespace optE { extern BooleanOptionKey const no_design_pdb_output; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
