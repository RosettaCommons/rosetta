// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/ms.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_ms_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_ms_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace ms { extern BooleanOptionKey const ms; }
namespace ms { extern BooleanOptionKey const share_data; }
namespace ms { extern BooleanOptionKey const verbose; }
namespace ms { extern BooleanOptionKey const restrict_to_canonical; }
namespace ms { extern IntegerOptionKey const pop_from_ss; }
namespace ms { extern IntegerOptionKey const pop_size; }
namespace ms { extern IntegerOptionKey const generations; }
namespace ms { extern IntegerOptionKey const num_packs; }
namespace ms { extern IntegerOptionKey const numresults; }
namespace ms { extern RealOptionKey const anchor_offset; }
namespace ms { extern RealOptionKey const Boltz_temp; }
namespace ms { extern RealOptionKey const mutate_rate; }
namespace ms { extern RealOptionKey const fraction_by_recombination; }
namespace ms { namespace checkpoint { extern BooleanOptionKey const checkpoint; } }
namespace ms { namespace checkpoint { extern StringOptionKey const prefix; } }
namespace ms { namespace checkpoint { extern IntegerOptionKey const interval; } }
namespace ms { namespace checkpoint { extern BooleanOptionKey const gz; } }
namespace ms { namespace checkpoint { extern BooleanOptionKey const rename; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
