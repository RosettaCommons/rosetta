// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/willmatch.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_willmatch_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_willmatch_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace willmatch { extern BooleanOptionKey const willmatch; }
namespace willmatch { extern RealOptionKey const arg_dun_th; }
namespace willmatch { extern RealOptionKey const asp_dun_th; }
namespace willmatch { extern RealOptionKey const glu_dun_th; }
namespace willmatch { extern RealOptionKey const lys_dun_th; }
namespace willmatch { extern BooleanOptionKey const usecache; }
namespace willmatch { extern StringVectorOptionKey const write_reduced_matchset; }
namespace willmatch { extern RealOptionKey const interface_size; }
namespace willmatch { extern RealOptionKey const max_dis_any; }
namespace willmatch { extern RealOptionKey const max_dis_all; }
namespace willmatch { extern RealOptionKey const max_dis_hb; }
namespace willmatch { extern RealOptionKey const min_dis_hb; }
namespace willmatch { extern RealOptionKey const max_dis_hb_colinear; }
namespace willmatch { extern RealOptionKey const max_dis_metal; }
namespace willmatch { extern RealOptionKey const max_ang_metal; }
namespace willmatch { extern RealOptionKey const clash_dis; }
namespace willmatch { extern RealOptionKey const c2_linker_dist; }
namespace willmatch { extern RealOptionKey const identical_match_dis; }
namespace willmatch { extern RealOptionKey const chi1_increment; }
namespace willmatch { extern RealOptionKey const chi2_increment; }
namespace willmatch { extern RealOptionKey const c2_symm_increment; }
namespace willmatch { extern RealOptionKey const cb_sasa_thresh; }
namespace willmatch { extern BooleanOptionKey const design_interface; }
namespace willmatch { extern FileOptionKey const chilist; }
namespace willmatch { extern FileOptionKey const fixed_res; }
namespace willmatch { extern FileOptionKey const native1; }
namespace willmatch { extern FileOptionKey const native2; }
namespace willmatch { extern FileOptionKey const exclude_res1; }
namespace willmatch { extern FileOptionKey const exclude_res2; }
namespace willmatch { extern FileOptionKey const taglist; }
namespace willmatch { extern IntegerVectorOptionKey const residues; }
namespace willmatch { extern BooleanOptionKey const symmetry_d2; }
namespace willmatch { extern BooleanOptionKey const symmetry_c2_dock; }
namespace willmatch { extern IntegerVectorOptionKey const splitwork; }
namespace willmatch { extern BooleanOptionKey const exclude_ala; }
namespace willmatch { extern RealOptionKey const match_overlap_dis; }
namespace willmatch { extern RealOptionKey const match_overlap_ang; }
namespace willmatch { extern IntegerVectorOptionKey const forbid_residues; }
namespace willmatch { extern RealVectorOptionKey const poi; }
namespace willmatch { extern RealOptionKey const poidis; }
namespace willmatch { extern BooleanOptionKey const homodimer; }
namespace willmatch { extern RealOptionKey const fa_dun_thresh; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
