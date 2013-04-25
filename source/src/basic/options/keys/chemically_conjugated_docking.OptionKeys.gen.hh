// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/chemically_conjugated_docking.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_chemically_conjugated_docking_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_chemically_conjugated_docking_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace chemically_conjugated_docking { extern BooleanOptionKey const chemically_conjugated_docking; }
namespace chemically_conjugated_docking { extern FileOptionKey const UBQpdb; }
namespace chemically_conjugated_docking { extern FileOptionKey const E2pdb; }
namespace chemically_conjugated_docking { extern IntegerOptionKey const E2_residue; }
namespace chemically_conjugated_docking { extern RealOptionKey const SASAfilter; }
namespace chemically_conjugated_docking { extern RealOptionKey const scorefilter; }
namespace chemically_conjugated_docking { extern BooleanOptionKey const publication; }
namespace chemically_conjugated_docking { extern IntegerOptionKey const n_tail_res; }
namespace chemically_conjugated_docking { extern BooleanOptionKey const two_ubiquitins; }
namespace chemically_conjugated_docking { extern FileVectorOptionKey const extra_bodies; }
namespace chemically_conjugated_docking { extern IntegerOptionKey const UBQ2_lys; }
namespace chemically_conjugated_docking { extern FileOptionKey const UBQ2_pdb; }
namespace chemically_conjugated_docking { extern BooleanOptionKey const dont_minimize_omega; }
namespace chemically_conjugated_docking { extern BooleanOptionKey const pdz; }
namespace chemically_conjugated_docking { extern FileOptionKey const GTPasepdb; }
namespace chemically_conjugated_docking { extern IntegerOptionKey const GTPase_residue; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
