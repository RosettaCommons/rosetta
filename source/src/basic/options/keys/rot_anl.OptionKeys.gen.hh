// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/rot_anl.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_rot_anl_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_rot_anl_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace rot_anl { extern BooleanOptionKey const rot_anl; }
namespace rot_anl { extern StringOptionKey const tag; }
namespace rot_anl { extern BooleanOptionKey const premin; }
namespace rot_anl { extern BooleanOptionKey const min; }
namespace rot_anl { extern BooleanOptionKey const diff_to_min; }
namespace rot_anl { extern BooleanOptionKey const repack; }
namespace rot_anl { extern BooleanOptionKey const rtmin; }
namespace rot_anl { extern BooleanOptionKey const scmove; }
namespace rot_anl { extern BooleanOptionKey const design; }
namespace rot_anl { extern RealOptionKey const score_tol; }
namespace rot_anl { extern RealOptionKey const rmsd_tol; }
namespace rot_anl { extern BooleanOptionKey const dump_pdb; }
namespace rot_anl { extern IntegerOptionKey const nloop_scmove; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
