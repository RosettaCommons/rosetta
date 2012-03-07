// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/robert.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_robert_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_robert_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace robert { extern BooleanOptionKey const robert; }
namespace robert { extern StringOptionKey const pairdata_input_pdb_list; }
namespace robert { extern RealOptionKey const pcs_maxsub_filter; }
namespace robert { extern RealOptionKey const pcs_maxsub_rmsd; }
namespace robert { extern BooleanOptionKey const pcs_dump_cluster; }
namespace robert { extern RealOptionKey const pcs_cluster_coverage; }
namespace robert { extern BooleanOptionKey const pcs_cluster_lowscoring; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
