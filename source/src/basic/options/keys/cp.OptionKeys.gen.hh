// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/cp.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_cp_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_cp_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace cp { extern BooleanOptionKey const cp; }
namespace cp { extern RealOptionKey const cutoff; }
namespace cp { extern StringOptionKey const minimizer; }
namespace cp { extern StringOptionKey const relax_sfxn; }
namespace cp { extern StringOptionKey const pack_sfxn; }
namespace cp { extern RealOptionKey const minimizer_tol; }
namespace cp { extern StringOptionKey const minimizer_score_fxn; }
namespace cp { extern StringOptionKey const output; }
namespace cp { extern IntegerOptionKey const ncycles; }
namespace cp { extern IntegerOptionKey const max_failures; }
namespace cp { extern BooleanOptionKey const print_reports; }
namespace cp { extern StringOptionKey const vipReportFile; }
namespace cp { extern StringOptionKey const exclude_file; }
namespace cp { extern StringOptionKey const relax_mover; }
namespace cp { extern BooleanOptionKey const skip_relax; }
namespace cp { extern BooleanOptionKey const local_relax; }
namespace cp { extern BooleanOptionKey const print_intermediate_pdbs; }
namespace cp { extern BooleanOptionKey const use_unrelaxed_starting_points; }
namespace cp { extern BooleanOptionKey const easy_vip_acceptance; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
