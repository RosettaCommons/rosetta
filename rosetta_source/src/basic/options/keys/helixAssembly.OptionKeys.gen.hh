// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/helixAssembly.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_helixAssembly_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_helixAssembly_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace helixAssembly { extern BooleanOptionKey const helixAssembly; }
namespace helixAssembly { extern FileOptionKey const query_structure_path; }
namespace helixAssembly { extern IntegerOptionKey const frag1_start; }
namespace helixAssembly { extern IntegerOptionKey const frag1_end; }
namespace helixAssembly { extern IntegerOptionKey const frag2_start; }
namespace helixAssembly { extern IntegerOptionKey const frag2_end; }
namespace helixAssembly { extern IntegerOptionKey const minimum_helix_contacts; }
namespace helixAssembly { extern IntegerOptionKey const helices_to_add; }
namespace helixAssembly { extern RealOptionKey const single_helix_rmsd_cutoff; }
namespace helixAssembly { extern RealOptionKey const helix_pair_rmsd_cutoff; }
namespace helixAssembly { extern RealOptionKey const helix_cap_dist_cutoff; }
namespace helixAssembly { extern RealOptionKey const helix_contact_dist_cutoff; }
namespace helixAssembly { extern IntegerOptionKey const min_helix_size; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
