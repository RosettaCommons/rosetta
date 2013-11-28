// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/edensity.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_edensity_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_edensity_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace edensity { extern BooleanOptionKey const edensity; }
namespace edensity { extern BooleanOptionKey const debug; }
namespace edensity { extern StringOptionKey const mapfile; }
namespace edensity { extern RealOptionKey const mapreso; }
namespace edensity { extern RealOptionKey const grid_spacing; }
namespace edensity { extern RealOptionKey const centroid_density_mass; }
namespace edensity { extern IntegerOptionKey const sliding_window; }
namespace edensity { extern BooleanOptionKey const cryoem_scatterers; }
namespace edensity { extern RealOptionKey const force_apix; }
namespace edensity { extern RealOptionKey const fastdens_wt; }
namespace edensity { extern RealVectorOptionKey const fastdens_params; }
namespace edensity { extern BooleanOptionKey const legacy_fastdens_score; }
namespace edensity { extern RealOptionKey const sliding_window_wt; }
namespace edensity { extern BooleanOptionKey const score_sliding_window_context; }
namespace edensity { extern RealOptionKey const whole_structure_ca_wt; }
namespace edensity { extern RealOptionKey const whole_structure_allatom_wt; }
namespace edensity { extern BooleanOptionKey const no_edens_in_minimizer; }
namespace edensity { extern BooleanOptionKey const debug_derivatives; }
namespace edensity { extern StringOptionKey const realign; }
namespace edensity { extern StringOptionKey const membrane_axis; }
namespace edensity { extern RealOptionKey const atom_mask; }
namespace edensity { extern RealOptionKey const ca_mask; }
namespace edensity { extern BooleanOptionKey const score_symm_complex; }
namespace edensity { extern RealOptionKey const sc_scaling; }
namespace edensity { extern IntegerOptionKey const n_kbins; }
namespace edensity { extern RealOptionKey const render_sigma; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
