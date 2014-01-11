// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/pH.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_pH_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_pH_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace pH { extern BooleanOptionKey const pH; }
namespace pH { extern BooleanOptionKey const pH_mode; }
namespace pH { extern BooleanOptionKey const keep_input_protonation_state; }
namespace pH { extern RealOptionKey const value_pH; }
namespace pH { namespace calc_pka { extern BooleanOptionKey const calc_pka; } }
namespace pH { namespace calc_pka { extern BooleanOptionKey const pka_all; } }
namespace pH { namespace calc_pka { extern RealVectorOptionKey const pka_for_resnos; } }
namespace pH { namespace calc_pka { extern StringOptionKey const pka_for_chainno; } }
namespace pH { namespace calc_pka { extern BooleanOptionKey const pH_neighbor_pack; } }
namespace pH { namespace calc_pka { extern RealOptionKey const pka_rad; } }
namespace pH { namespace calc_pka { extern BooleanOptionKey const pH_prepack; } }
namespace pH { namespace calc_pka { extern BooleanOptionKey const pH_relax; } }
namespace pH { namespace calc_pka { extern BooleanOptionKey const rotamer_prot_stats; } }
namespace pH { extern FileVectorOptionKey const pH_unbound; }
namespace pH { extern BooleanOptionKey const output_raw_scores; }
namespace pH { extern BooleanOptionKey const pre_process; }
namespace pH { extern StringOptionKey const cognate_partners; }
namespace pH { extern FileOptionKey const cognate_pdb; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
