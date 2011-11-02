// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/remodel.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_remodel_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_remodel_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>


namespace basic {
namespace options {
namespace OptionKeys {

namespace remodel { extern BooleanOptionKey const remodel; }
namespace remodel { extern BooleanOptionKey const help; }
namespace remodel { extern BooleanOptionKey const autopilot; }
namespace remodel { extern FileOptionKey const blueprint; }
namespace remodel { extern FileOptionKey const cstfile; }
namespace remodel { extern IntegerOptionKey const num_trajectory; }
namespace remodel { extern IntegerOptionKey const save_top; }
namespace remodel { extern BooleanOptionKey const swap_refine_confirm_protocols; }
namespace remodel { extern IntegerOptionKey const num_frag_moves; }
namespace remodel { extern BooleanOptionKey const bypass_fragments; }
namespace remodel { extern BooleanOptionKey const use_same_length_fragments; }
namespace remodel { extern BooleanOptionKey const enable_ligand_aa; }
namespace remodel { extern BooleanOptionKey const no_jumps; }
namespace remodel { extern BooleanOptionKey const backrub; }
namespace remodel { extern BooleanOptionKey const use_blueprint_sequence ; }
namespace remodel { extern BooleanOptionKey const randomize_equivalent_fragments ; }
namespace remodel { extern BooleanOptionKey const quick_and_dirty ; }
namespace remodel { extern BooleanOptionKey const checkpoint ; }
namespace remodel { extern BooleanOptionKey const use_ccd_refine ; }
namespace remodel { extern BooleanOptionKey const use_pose_relax ; }
namespace remodel { extern BooleanOptionKey const use_dssp_assignment; }
namespace remodel { extern BooleanOptionKey const keep_jumps_in_minimizer ; }
namespace remodel { extern FileOptionKey const output_fragfiles; }
namespace remodel { extern FileOptionKey const read_fragfile; }
namespace remodel { extern StringOptionKey const generic_aa; }
namespace remodel { extern RealOptionKey const cluster_radius; }
namespace remodel { extern BooleanOptionKey const use_clusters; }
namespace remodel { extern BooleanOptionKey const run_confirmation; }
namespace remodel { extern BooleanOptionKey const cluster_on_entire_pose; }
namespace remodel { extern IntegerOptionKey const collect_clustered_top; }
namespace remodel { extern IntegerOptionKey const dr_cycles; }
namespace remodel { extern IntegerOptionKey const two_chain_tree; }
namespace remodel { extern IntegerOptionKey const repeat_structuer; }
namespace remodel { namespace domainFusion { extern BooleanOptionKey const domainFusion; } }
namespace remodel { namespace domainFusion { extern FileOptionKey const insert_segment_from_pdb; } }
namespace remodel { extern RealOptionKey const vdw; }
namespace remodel { extern RealOptionKey const rama; }
namespace remodel { extern RealOptionKey const cbeta; }
namespace remodel { extern RealOptionKey const cenpack; }
namespace remodel { extern RealOptionKey const hb_lrbb; }
namespace remodel { extern RealOptionKey const hb_srbb; }
namespace remodel { extern RealOptionKey const rg; }
namespace remodel { extern RealOptionKey const rsigma; }
namespace remodel { extern RealOptionKey const ss_pair; }
namespace remodel { extern BooleanOptionKey const build_disulf; }
namespace remodel { extern IntegerOptionKey const max_disulf_allowed; }
namespace remodel { extern RealOptionKey const match_rt_limit; }
namespace remodel { extern IntegerVectorOptionKey const disulf_landing_range; }
namespace remodel { namespace design { extern BooleanOptionKey const design; } }
namespace remodel { namespace design { extern BooleanOptionKey const no_design ; } }
namespace remodel { namespace design { extern BooleanOptionKey const silent; } }
namespace remodel { namespace design { extern BooleanOptionKey const skip_partial; } }
namespace remodel { namespace design { extern BooleanOptionKey const design_neighbors; } }
namespace remodel { namespace design { extern BooleanOptionKey const find_neighbors; } }
namespace remodel { extern BooleanOptionKey const rank_by_bsasa; }
namespace remodel { namespace RemodelLoopMover { extern BooleanOptionKey const RemodelLoopMover; } }
namespace remodel { namespace RemodelLoopMover { extern IntegerOptionKey const max_linear_chainbreak; } }
namespace remodel { namespace RemodelLoopMover { extern BooleanOptionKey const randomize_loops; } }
namespace remodel { namespace RemodelLoopMover { extern IntegerOptionKey const allowed_closure_attempts; } }
namespace remodel { namespace RemodelLoopMover { extern IntegerOptionKey const simultaneous_cycles; } }
namespace remodel { namespace RemodelLoopMover { extern IntegerOptionKey const independent_cycles; } }
namespace remodel { namespace RemodelLoopMover { extern IntegerOptionKey const boost_closure_cycles; } }
namespace remodel { namespace RemodelLoopMover { extern BooleanOptionKey const force_cutting_N; } }
namespace remodel { namespace RemodelLoopMover { extern BooleanOptionKey const bypass_closure; } }
namespace remodel { namespace RemodelLoopMover { extern RealOptionKey const temperature; } }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
