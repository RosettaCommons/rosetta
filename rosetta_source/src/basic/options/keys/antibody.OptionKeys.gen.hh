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
namespace antibody { extern BooleanOptionKey const model_h3; }
namespace antibody { extern BooleanOptionKey const snugfit; }
namespace antibody { extern BooleanOptionKey const refine_h3; }
namespace antibody { extern BooleanOptionKey const h3_filter; }
namespace antibody { extern BooleanOptionKey const cter_insert; }
namespace antibody { extern BooleanOptionKey const flank_residue_min; }
namespace antibody { extern BooleanOptionKey const sc_min; }
namespace antibody { extern BooleanOptionKey const rt_min; }
namespace antibody { extern StringOptionKey const remodel; }
namespace antibody { extern StringOptionKey const refine; }
namespace antibody { extern BooleanOptionKey const snug_loops; }
namespace antibody { extern FileOptionKey const input_fv; }
namespace antibody { extern BooleanOptionKey const camelid; }
namespace antibody { extern BooleanOptionKey const camelid_constraints; }
namespace antibody { extern StringOptionKey const numbering_scheme; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
