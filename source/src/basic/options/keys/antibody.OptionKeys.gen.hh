// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/antibody.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_antibody_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_antibody_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace antibody { extern BooleanOptionKey const antibody; }
namespace antibody { extern BooleanOptionKey const graft_l1; }
namespace antibody { extern StringOptionKey const l1_template; }
namespace antibody { extern BooleanOptionKey const graft_l2; }
namespace antibody { extern StringOptionKey const l2_template; }
namespace antibody { extern BooleanOptionKey const graft_l3; }
namespace antibody { extern StringOptionKey const l3_template; }
namespace antibody { extern BooleanOptionKey const graft_h1; }
namespace antibody { extern StringOptionKey const h1_template; }
namespace antibody { extern BooleanOptionKey const graft_h2; }
namespace antibody { extern StringOptionKey const h2_template; }
namespace antibody { extern BooleanOptionKey const graft_h3; }
namespace antibody { extern StringOptionKey const h3_template; }
namespace antibody { extern BooleanOptionKey const h3_no_stem_graft; }
namespace antibody { extern BooleanOptionKey const packonly_after_graft; }
namespace antibody { extern BooleanOptionKey const stem_optimize; }
namespace antibody { extern IntegerOptionKey const stem_optimize_size; }
namespace antibody { extern StringOptionKey const preprocessing_script_version; }
namespace antibody { extern BooleanOptionKey const model_h3; }
namespace antibody { extern BooleanOptionKey const snugfit; }
namespace antibody { extern BooleanOptionKey const refine_h3; }
namespace antibody { extern BooleanOptionKey const h3_filter; }
namespace antibody { extern RealOptionKey const h3_filter_tolerance; }
namespace antibody { extern BooleanOptionKey const cter_insert; }
namespace antibody { extern BooleanOptionKey const flank_residue_min; }
namespace antibody { extern BooleanOptionKey const sc_min; }
namespace antibody { extern BooleanOptionKey const rt_min; }
namespace antibody { extern BooleanOptionKey const bad_nter; }
namespace antibody { extern BooleanOptionKey const extend_h3_before_modeling; }
namespace antibody { extern BooleanOptionKey const idealize_h3_stems_before_modeling; }
namespace antibody { extern StringOptionKey const remodel; }
namespace antibody { extern StringOptionKey const refine; }
namespace antibody { extern StringOptionKey const centroid_refine; }
namespace antibody { extern BooleanOptionKey const constrain_cter; }
namespace antibody { extern BooleanOptionKey const constrain_vlvh_qq; }
namespace antibody { extern BooleanOptionKey const snug_loops; }
namespace antibody { extern FileOptionKey const input_fv; }
namespace antibody { extern BooleanOptionKey const camelid; }
namespace antibody { extern BooleanOptionKey const camelid_constraints; }
namespace antibody { extern StringOptionKey const numbering_scheme; }
namespace antibody { namespace design { extern BooleanOptionKey const design; } }
namespace antibody { namespace design { extern StringOptionKey const instructions; } }
namespace antibody { namespace design { extern StringOptionKey const antibody_database; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_graft_design; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_post_graft_design_modeling; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_sequence_design; } }
namespace antibody { namespace design { extern BooleanOptionKey const do_post_design_modeling; } }
namespace antibody { namespace design { extern IntegerOptionKey const graft_rounds; } }
namespace antibody { namespace design { extern IntegerOptionKey const top_graft_designs; } }
namespace antibody { namespace design { extern BooleanOptionKey const initial_perturb; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_deterministic; } }
namespace antibody { namespace design { extern BooleanOptionKey const dump_post_graft_designs; } }
namespace antibody { namespace design { extern RealOptionKey const interface_dis; } }
namespace antibody { namespace design { extern RealOptionKey const neighbor_dis; } }
namespace antibody { namespace design { extern BooleanOptionKey const dock_post_graft; } }
namespace antibody { namespace design { extern BooleanOptionKey const pack_post_graft; } }
namespace antibody { namespace design { extern BooleanOptionKey const rb_min_post_graft; } }
namespace antibody { namespace design { extern BooleanOptionKey const design_post_graft; } }
namespace antibody { namespace design { extern IntegerOptionKey const dock_rounds; } }
namespace antibody { namespace design { extern StringOptionKey const ab_dock_chains; } }
namespace antibody { namespace design { extern StringOptionKey const design_method; } }
namespace antibody { namespace design { extern IntegerOptionKey const design_rounds; } }
namespace antibody { namespace design { extern StringOptionKey const design_scorefxn; } }
namespace antibody { namespace design { extern BooleanOptionKey const benchmark_basic_design; } }
namespace antibody { namespace design { extern BooleanOptionKey const use_filters; } }
namespace antibody { namespace design { extern IntegerOptionKey const stats_cutoff; } }
namespace antibody { namespace design { extern IntegerOptionKey const sample_zero_probs_at; } }
namespace antibody { namespace design { extern BooleanOptionKey const conservative_h3_design; } }
namespace antibody { namespace design { extern BooleanOptionKey const turn_conservation; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
