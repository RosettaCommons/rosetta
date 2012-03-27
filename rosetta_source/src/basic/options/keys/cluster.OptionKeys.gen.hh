// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/cluster.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_cluster_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_cluster_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace cluster { extern BooleanOptionKey const cluster; }
namespace cluster { extern BooleanOptionKey const lite; }
namespace cluster { extern RealOptionKey const input_score_filter; }
namespace cluster { extern RealOptionKey const output_score_filter; }
namespace cluster { extern IntegerVectorOptionKey const exclude_res; }
namespace cluster { extern RealOptionKey const thinout_factor; }
namespace cluster { extern IntegerOptionKey const max_cluster_seeds; }
namespace cluster { extern RealOptionKey const radius; }
namespace cluster { extern IntegerOptionKey const limit_cluster_size; }
namespace cluster { extern RealOptionKey const limit_cluster_size_percent; }
namespace cluster { extern RealOptionKey const random_limit_cluster_size_percent; }
namespace cluster { extern IntegerOptionKey const limit_clusters; }
namespace cluster { extern IntegerOptionKey const limit_total_structures; }
namespace cluster { extern IntegerOptionKey const max_total_cluster; }
namespace cluster { extern BooleanOptionKey const gdtmm; }
namespace cluster { extern BooleanOptionKey const sort_groups_by_energy; }
namespace cluster { extern BooleanOptionKey const sort_groups_by_size; }
namespace cluster { extern BooleanOptionKey const remove_singletons; }
namespace cluster { extern BooleanOptionKey const export_only_low; }
namespace cluster { extern BooleanOptionKey const remove_highest_energy_member; }
namespace cluster { extern BooleanOptionKey const idealize_final_structures; }
namespace cluster { extern IntegerOptionKey const limit_dist_matrix; }
namespace cluster { extern BooleanOptionKey const make_ensemble_cst; }
namespace cluster { extern BooleanOptionKey const hotspot_hash; }
namespace cluster { extern BooleanOptionKey const loops; }
namespace cluster { extern RealOptionKey const population_weight; }
namespace cluster { extern StringOptionKey const template_scores; }
namespace cluster { extern IntegerOptionKey const K_level; }
namespace cluster { extern RealVectorOptionKey const K_radius; }
namespace cluster { extern IntegerVectorOptionKey const K_n_cluster; }
namespace cluster { extern StringVectorOptionKey const K_style; }
namespace cluster { extern RealOptionKey const K_threshold; }
namespace cluster { extern IntegerOptionKey const K_n_sub; }
namespace cluster { extern IntegerOptionKey const K_deque_size; }
namespace cluster { extern IntegerOptionKey const K_deque_level; }
namespace cluster { extern BooleanOptionKey const K_redundant; }
namespace cluster { extern BooleanOptionKey const K_not_fit_xyz; }
namespace cluster { extern BooleanOptionKey const K_save_headers; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
