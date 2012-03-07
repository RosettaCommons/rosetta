// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/flexPepDocking.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_flexPepDocking_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_flexPepDocking_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace flexPepDocking { extern BooleanOptionKey const flexPepDocking; }
namespace flexPepDocking { extern StringOptionKey const params_file; }
namespace flexPepDocking { extern IntegerOptionKey const peptide_anchor; }
namespace flexPepDocking { extern StringOptionKey const receptor_chain; }
namespace flexPepDocking { extern StringOptionKey const peptide_chain; }
namespace flexPepDocking { extern BooleanOptionKey const pep_fold_only; }
namespace flexPepDocking { extern BooleanOptionKey const lowres_abinitio; }
namespace flexPepDocking { extern BooleanOptionKey const lowres_preoptimize; }
namespace flexPepDocking { extern BooleanOptionKey const flexPepDockingMinimizeOnly; }
namespace flexPepDocking { extern BooleanOptionKey const extend_peptide; }
namespace flexPepDocking { extern BooleanOptionKey const pep_refine; }
namespace flexPepDocking { extern BooleanOptionKey const rbMCM; }
namespace flexPepDocking { extern BooleanOptionKey const torsionsMCM; }
namespace flexPepDocking { extern BooleanOptionKey const peptide_loop_model; }
namespace flexPepDocking { extern BooleanOptionKey const backrub_peptide; }
namespace flexPepDocking { extern BooleanOptionKey const boost_fa_atr; }
namespace flexPepDocking { extern BooleanOptionKey const ramp_fa_rep; }
namespace flexPepDocking { extern BooleanOptionKey const ramp_rama; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_score_only; }
namespace flexPepDocking { extern FileOptionKey const ref_startstruct; }
namespace flexPepDocking { extern BooleanOptionKey const use_cen_score; }
namespace flexPepDocking { extern BooleanOptionKey const design_peptide; }
namespace flexPepDocking { extern IntegerOptionKey const rep_ramp_cycles; }
namespace flexPepDocking { extern IntegerOptionKey const mcm_cycles; }
namespace flexPepDocking { extern RealOptionKey const random_phi_psi_preturbation; }
namespace flexPepDocking { extern RealOptionKey const smove_angle_range; }
namespace flexPepDocking { extern BooleanOptionKey const min_receptor_bb; }
namespace flexPepDocking { extern RealOptionKey const random_trans_start; }
namespace flexPepDocking { extern RealOptionKey const random_rot_start; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_prepack; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_noprepack1; }
namespace flexPepDocking { extern BooleanOptionKey const flexpep_noprepack2; }
namespace flexPepDocking { extern RealOptionKey const score_filter; }
namespace flexPepDocking { extern IntegerOptionKey const hb_filter; }
namespace flexPepDocking { extern IntegerOptionKey const hotspot_filter; }
namespace flexPepDocking { extern StringOptionKey const frag5; }
namespace flexPepDocking { extern RealOptionKey const frag9_weight; }
namespace flexPepDocking { extern RealOptionKey const frag5_weight; }
namespace flexPepDocking { extern RealOptionKey const frag3_weight; }
namespace flexPepDocking { extern BooleanOptionKey const pSer2Asp_centroid; }
namespace flexPepDocking { extern BooleanOptionKey const pSer2Glu_centroid; }
namespace flexPepDocking { extern BooleanOptionKey const dumpPDB_abinitio; }
namespace flexPepDocking { extern BooleanOptionKey const dumpPDB_lowres; }
namespace flexPepDocking { extern BooleanOptionKey const dumpPDB_hires; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
