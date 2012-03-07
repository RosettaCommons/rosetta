// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/hotspot.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_hotspot_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_hotspot_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace hotspot { extern BooleanOptionKey const hotspot; }
namespace hotspot { extern BooleanOptionKey const allow_gly; }
namespace hotspot { extern BooleanOptionKey const allow_proline; }
namespace hotspot { extern BooleanOptionKey const benchmark; }
namespace hotspot { extern StringVectorOptionKey const residue; }
namespace hotspot { extern FileOptionKey const hashfile; }
namespace hotspot { extern FileOptionKey const target; }
namespace hotspot { extern IntegerOptionKey const target_res; }
namespace hotspot { extern RealOptionKey const target_dist; }
namespace hotspot { extern FileOptionKey const density; }
namespace hotspot { extern FileOptionKey const weighted_density; }
namespace hotspot { extern FileOptionKey const rms_target; }
namespace hotspot { extern FileOptionKey const rms_hotspot; }
namespace hotspot { extern IntegerOptionKey const rms_hotspot_res; }
namespace hotspot { extern BooleanOptionKey const rescore; }
namespace hotspot { extern RealOptionKey const threshold; }
namespace hotspot { extern BooleanOptionKey const sc_only; }
namespace hotspot { extern BooleanOptionKey const fxnal_group; }
namespace hotspot { extern BooleanOptionKey const cluster; }
namespace hotspot { extern BooleanOptionKey const colonyE; }
namespace hotspot { extern IntegerOptionKey const length; }
namespace hotspot { extern BooleanOptionKey const envhb; }
namespace hotspot { extern RealOptionKey const angle; }
namespace hotspot { extern IntegerOptionKey const angle_res; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
