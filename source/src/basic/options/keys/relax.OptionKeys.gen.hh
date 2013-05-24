// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/relax.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_relax_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_relax_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace relax { extern BooleanOptionKey const relax; }
namespace relax { extern BooleanOptionKey const fast; }
namespace relax { extern BooleanOptionKey const thorough; }
namespace relax { extern BooleanOptionKey const membrane; }
namespace relax { extern BooleanOptionKey const centroid_mode; }
namespace relax { extern IntegerOptionKey const default_repeats; }
namespace relax { extern BooleanOptionKey const dualspace; }
namespace relax { extern BooleanOptionKey const ramady; }
namespace relax { extern RealOptionKey const ramady_rms_limit; }
namespace relax { extern RealOptionKey const ramady_cutoff; }
namespace relax { extern IntegerOptionKey const ramady_max_rebuild; }
namespace relax { extern BooleanOptionKey const ramady_force; }
namespace relax { extern FileOptionKey const script; }
namespace relax { extern IntegerOptionKey const script_max_accept; }
namespace relax { extern BooleanOptionKey const superimpose_to_native; }
namespace relax { extern FileOptionKey const superimpose_to_file; }
namespace relax { extern BooleanOptionKey const constrain_relax_to_native_coords; }
namespace relax { extern BooleanOptionKey const constrain_relax_to_start_coords; }
namespace relax { extern BooleanOptionKey const coord_constrain_sidechains; }
namespace relax { extern RealOptionKey const sc_cst_maxdist; }
namespace relax { extern BooleanOptionKey const limit_aroma_chi2; }
namespace relax { extern BooleanOptionKey const respect_resfile; }
namespace relax { extern BooleanOptionKey const bb_move; }
namespace relax { extern BooleanOptionKey const chi_move; }
namespace relax { extern BooleanOptionKey const jump_move; }
namespace relax { extern BooleanOptionKey const dna_move; }
namespace relax { extern BooleanOptionKey const fix_omega; }
namespace relax { extern BooleanOptionKey const minimize_bond_lengths; }
namespace relax { extern BooleanOptionKey const minimize_bond_angles; }
namespace relax { extern IntegerOptionKey const minimize_bondlength_subset; }
namespace relax { extern IntegerOptionKey const minimize_bondangle_subset; }
namespace relax { extern StringOptionKey const min_type; }
namespace relax { extern BooleanOptionKey const cartesian; }
namespace relax { extern RealOptionKey const chainbreak_weight; }
namespace relax { extern RealOptionKey const linear_chainbreak_weight; }
namespace relax { extern RealOptionKey const overlap_chainbreak_weight; }
namespace relax { extern BooleanOptionKey const classic; }
namespace relax { extern FileOptionKey const sequence_file; }
namespace relax { extern BooleanOptionKey const quick; }
namespace relax { extern BooleanOptionKey const sequence; }
namespace relax { extern IntegerOptionKey const minirelax_repeats; }
namespace relax { extern RealOptionKey const minirelax_sdev; }
namespace relax { extern BooleanOptionKey const wobblemoves; }
namespace relax { extern FileOptionKey const constrain_relax_segments; }
namespace relax { extern RealOptionKey const coord_cst_width; }
namespace relax { extern RealOptionKey const coord_cst_stdev; }
namespace relax { extern BooleanOptionKey const ramp_constraints; }
namespace relax { extern RealOptionKey const energycut; }
namespace relax { extern BooleanOptionKey const mini; }
namespace relax { extern IntegerOptionKey const stage1_ramp_cycles; }
namespace relax { extern IntegerOptionKey const stage1_ramp_inner_cycles; }
namespace relax { extern IntegerOptionKey const stage2_repack_period; }
namespace relax { extern IntegerOptionKey const stage2_cycles; }
namespace relax { extern RealOptionKey const min_tolerance; }
namespace relax { extern IntegerOptionKey const stage3_cycles; }
namespace relax { extern RealOptionKey const cycle_ratio; }
namespace relax { extern RealOptionKey const filter_stage2_beginning; }
namespace relax { extern RealOptionKey const filter_stage2_quarter; }
namespace relax { extern RealOptionKey const filter_stage2_half; }
namespace relax { extern RealOptionKey const filter_stage2_end; }
namespace relax { namespace centroid { extern BooleanOptionKey const centroid; } }
namespace relax { namespace centroid { extern StringOptionKey const weights; } }
namespace relax { namespace centroid { extern BooleanOptionKey const ramp_vdw; } }
namespace relax { namespace centroid { extern BooleanOptionKey const ramp_rama; } }
namespace relax { namespace centroid { extern StringOptionKey const parameters; } }
namespace relax { namespace centroid { extern BooleanOptionKey const do_final_repack; } }
namespace relax { namespace centroid { extern BooleanOptionKey const increase_vdw_radii; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
