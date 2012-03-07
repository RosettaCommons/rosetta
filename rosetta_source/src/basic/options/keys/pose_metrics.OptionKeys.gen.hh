// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/pose_metrics.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_pose_metrics_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_pose_metrics_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace pose_metrics { extern BooleanOptionKey const pose_metrics; }
namespace pose_metrics { extern RealOptionKey const atomic_burial_cutoff; }
namespace pose_metrics { extern RealOptionKey const sasa_calculator_probe_radius; }
namespace pose_metrics { extern RealOptionKey const interface_cutoff; }
namespace pose_metrics { extern IntegerOptionKey const min_sequence_separation; }
namespace pose_metrics { extern RealOptionKey const contact_cutoffE; }
namespace pose_metrics { extern RealOptionKey const neighbor_by_distance_cutoff; }
namespace pose_metrics { extern RealOptionKey const inter_group_neighbors_cutoff; }
namespace pose_metrics { extern RealOptionKey const semiex_water_burial_cutoff; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
