// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/pocket_grid.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_pocket_grid_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_pocket_grid_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace pocket_grid { extern BooleanOptionKey const pocket_grid; }
namespace pocket_grid { extern RealOptionKey const pocket_grid_size; }
namespace pocket_grid { extern RealOptionKey const pocket_grid_size_x; }
namespace pocket_grid { extern RealOptionKey const pocket_grid_size_y; }
namespace pocket_grid { extern RealOptionKey const pocket_grid_size_z; }
namespace pocket_grid { extern RealOptionKey const pocket_grid_spacing; }
namespace pocket_grid { extern RealOptionKey const pocket_max_spacing; }
namespace pocket_grid { extern RealOptionKey const pocket_min_size; }
namespace pocket_grid { extern RealOptionKey const pocket_max_size; }
namespace pocket_grid { extern RealOptionKey const pocket_probe_radius; }
namespace pocket_grid { extern StringOptionKey const central_relax_pdb_num; }
namespace pocket_grid { extern IntegerOptionKey const pocket_ntrials; }
namespace pocket_grid { extern IntegerOptionKey const pocket_num_angles; }
namespace pocket_grid { extern BooleanOptionKey const pocket_side; }
namespace pocket_grid { extern BooleanOptionKey const pocket_dump_pdbs; }
namespace pocket_grid { extern BooleanOptionKey const pocket_dump_exemplars; }
namespace pocket_grid { extern BooleanOptionKey const pocket_filter_by_exemplar; }
namespace pocket_grid { extern BooleanOptionKey const pocket_dump_rama; }
namespace pocket_grid { extern BooleanOptionKey const pocket_restrict_size; }
namespace pocket_grid { extern BooleanOptionKey const pocket_ignore_buried; }
namespace pocket_grid { extern BooleanOptionKey const pocket_only_buried; }
namespace pocket_grid { extern BooleanOptionKey const pocket_psp; }
namespace pocket_grid { extern BooleanOptionKey const pocket_sps; }
namespace pocket_grid { extern BooleanOptionKey const pocket_search13; }
namespace pocket_grid { extern RealOptionKey const pocket_surface_score; }
namespace pocket_grid { extern RealOptionKey const pocket_surface_dist; }
namespace pocket_grid { extern RealOptionKey const pocket_buried_score; }
namespace pocket_grid { extern RealOptionKey const pocket_buried_dist; }
namespace pocket_grid { extern RealOptionKey const pocket_exemplar_vdw_pen; }
namespace pocket_grid { extern BooleanOptionKey const pocket_debug_output; }
namespace pocket_grid { extern BooleanOptionKey const print_grid; }
namespace pocket_grid { extern BooleanOptionKey const extend_eggshell; }
namespace pocket_grid { extern RealOptionKey const extend_eggshell_dist; }
namespace pocket_grid { extern RealOptionKey const extra_eggshell_dist; }
namespace pocket_grid { extern RealOptionKey const eggshell_dist; }
namespace pocket_grid { extern BooleanOptionKey const reduce_rays; }
namespace pocket_grid { extern BooleanOptionKey const pocket_static_grid; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
