// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/james.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_james_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_james_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace james { extern BooleanOptionKey const james; }
namespace james { extern IntegerOptionKey const min_seqsep; }
namespace james { extern StringVectorOptionKey const atom_names; }
namespace james { extern RealVectorOptionKey const dist_thresholds; }
namespace james { extern RealVectorOptionKey const torsion_thresholds; }
namespace james { extern RealOptionKey const sog_cutoff; }
namespace james { extern BooleanOptionKey const shift_sog_func; }
namespace james { extern StringOptionKey const min_type; }
namespace james { extern RealOptionKey const min_tol; }
namespace james { extern BooleanOptionKey const debug; }
namespace james { extern RealOptionKey const real; }
namespace james { extern IntegerOptionKey const n_designs; }
namespace james { extern BooleanOptionKey const awesome_mode; }
namespace james { extern IntegerOptionKey const n_clusters; }
namespace james { extern BooleanOptionKey const thread_unaligned; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
