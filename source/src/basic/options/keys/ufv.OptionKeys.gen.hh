// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/ufv.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_ufv_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_ufv_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace ufv { extern BooleanOptionKey const ufv; }
namespace ufv { extern IntegerOptionKey const left; }
namespace ufv { extern IntegerOptionKey const right; }
namespace ufv { extern StringOptionKey const ss; }
namespace ufv { extern StringOptionKey const aa_during_build; }
namespace ufv { extern StringOptionKey const aa_during_design_refine; }
namespace ufv { extern BooleanOptionKey const keep_junction_torsions; }
namespace ufv { extern FileOptionKey const ufv_loops; }
namespace ufv { extern BooleanOptionKey const use_fullmer; }
namespace ufv { extern StringOptionKey const centroid_loop_mover; }
namespace ufv { extern BooleanOptionKey const no_neighborhood_design; }
namespace ufv { extern IntegerOptionKey const dr_cycles; }
namespace ufv { extern StringOptionKey const centroid_sfx; }
namespace ufv { extern StringOptionKey const centroid_sfx_patch; }
namespace ufv { extern StringOptionKey const fullatom_sfx; }
namespace ufv { extern StringOptionKey const fullatom_sfx_patch; }
namespace ufv { namespace insert { extern BooleanOptionKey const insert; } }
namespace ufv { namespace insert { extern FileOptionKey const insert_pdb; } }
namespace ufv { namespace insert { extern FileOptionKey const attached_pdb; } }
namespace ufv { namespace insert { extern StringOptionKey const connection_scheme; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
