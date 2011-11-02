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

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>


namespace basic {
namespace options {
namespace OptionKeys {

namespace optE { extern BooleanOptionKey const optE; }
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
namespace optE { extern BooleanOptionKey const optE_no_protein_hack_elec; }
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
namespace optE { namespace rescore { extern BooleanOptionKey const rescore; } }
namespace optE { namespace rescore { extern FileOptionKey const weights; } }
namespace optE { namespace rescore { extern IntegerOptionKey const context_round; } }
namespace optE { namespace rescore { extern FileOptionKey const outlog; } }
namespace optE { namespace rescore { extern BooleanOptionKey const measure_sequence_recovery; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
