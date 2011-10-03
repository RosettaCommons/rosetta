// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/rdc.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_rdc_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_rdc_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace rdc { extern BooleanOptionKey const rdc; }
namespace rdc { extern BooleanOptionKey const correct_NH_length; }
namespace rdc { extern BooleanOptionKey const reduced_couplings; }
namespace rdc { extern FileOptionKey const weights; }
namespace rdc { extern RealOptionKey const iterate_weights; }
namespace rdc { extern FileOptionKey const segment_file; }
namespace rdc { extern StringOptionKey const segment_scoring_mode; }
namespace rdc { extern RealOptionKey const total_weight; }
namespace rdc { extern RealOptionKey const tensor_weight; }
namespace rdc { extern FileOptionKey const print_rdc_values; }
namespace rdc { extern RealOptionKey const iterate_tol; }
namespace rdc { extern BooleanOptionKey const iterate_reset; }
namespace rdc { extern FileOptionKey const dump_weight_trajectory; }
namespace rdc { extern RealVectorOptionKey const fix_normAzz; }
namespace rdc { extern FileOptionKey const select_residues_file; }
namespace rdc { extern StringOptionKey const fit_method; }
namespace rdc { extern RealVectorOptionKey const fixDa; }
namespace rdc { extern RealVectorOptionKey const fixR; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
