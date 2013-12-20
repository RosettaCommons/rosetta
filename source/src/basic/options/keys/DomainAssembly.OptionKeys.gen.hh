// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/DomainAssembly.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_DomainAssembly_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_DomainAssembly_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace DomainAssembly { extern BooleanOptionKey const DomainAssembly; }
namespace DomainAssembly { extern BooleanOptionKey const da_setup; }
namespace DomainAssembly { extern FileOptionKey const da_setup_option_file; }
namespace DomainAssembly { extern FileOptionKey const da_setup_output_pdb; }
namespace DomainAssembly { extern FileOptionKey const da_linker_file; }
namespace DomainAssembly { extern FileOptionKey const da_require_buried; }
namespace DomainAssembly { extern FileOptionKey const da_start_pdb; }
namespace DomainAssembly { extern BooleanOptionKey const run_fullatom; }
namespace DomainAssembly { extern BooleanOptionKey const run_centroid; }
namespace DomainAssembly { extern BooleanOptionKey const run_centroid_abinitio; }
namespace DomainAssembly { extern IntegerOptionKey const da_nruns; }
namespace DomainAssembly { extern IntegerOptionKey const da_start_pdb_num; }
namespace DomainAssembly { extern FileOptionKey const da_linker_file_rna; }
namespace DomainAssembly { extern StringOptionKey const residues_repack_only; }
namespace DomainAssembly { extern FileOptionKey const da_eval_pose_map; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
