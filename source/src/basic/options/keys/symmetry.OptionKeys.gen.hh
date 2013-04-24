// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/symmetry.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_symmetry_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_symmetry_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace symmetry { extern BooleanOptionKey const symmetry; }
namespace symmetry { extern StringOptionKey const symmetry_definition; }
namespace symmetry { extern RealOptionKey const reweight_symm_interactions; }
namespace symmetry { extern BooleanOptionKey const initialize_rigid_body_dofs; }
namespace symmetry { extern BooleanOptionKey const detect_bonds; }
namespace symmetry { extern RealVectorOptionKey const perturb_rigid_body_dofs; }
namespace symmetry { extern BooleanOptionKey const symmetric_rmsd; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
