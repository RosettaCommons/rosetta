// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/sasa.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_sasa_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_sasa_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace sasa { extern BooleanOptionKey const sasa; }
namespace sasa { extern StringOptionKey const method; }
namespace sasa { extern BooleanOptionKey const include_hydrogens_explicitly; }
namespace sasa { extern RealOptionKey const probe_radius; }
namespace sasa { extern BooleanOptionKey const include_probe_radius_in_atom_radii; }
namespace sasa { extern BooleanOptionKey const include_only_C_S_in_hsasa; }
namespace sasa { extern BooleanOptionKey const exclude_polar_atoms_by_charge_in_hsasa; }
namespace sasa { extern RealOptionKey const polar_charge_cutoff; }
namespace sasa { extern StringOptionKey const implicit_hydrogen_radii_set; }
namespace sasa { extern StringOptionKey const explicit_hydrogen_radii_set; }
namespace sasa { extern BooleanOptionKey const use_legacy_behavior; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
