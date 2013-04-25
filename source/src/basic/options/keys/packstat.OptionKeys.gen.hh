// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/packstat.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_packstat_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_packstat_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace packstat { extern BooleanOptionKey const packstat; }
namespace packstat { extern BooleanOptionKey const include_water; }
namespace packstat { extern IntegerOptionKey const oversample; }
namespace packstat { extern BooleanOptionKey const packstat_pdb; }
namespace packstat { extern BooleanOptionKey const surface_accessibility; }
namespace packstat { extern BooleanOptionKey const residue_scores; }
namespace packstat { extern RealOptionKey const cavity_burial_probe_radius; }
namespace packstat { extern BooleanOptionKey const raw_stats; }
namespace packstat { extern IntegerOptionKey const threads; }
namespace packstat { extern RealOptionKey const cluster_min_volume; }
namespace packstat { extern RealOptionKey const min_surface_accessibility; }
namespace packstat { extern RealOptionKey const min_cluster_overlap; }
namespace packstat { extern RealOptionKey const min_cav_ball_radius; }
namespace packstat { extern RealOptionKey const max_cav_ball_radius; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
