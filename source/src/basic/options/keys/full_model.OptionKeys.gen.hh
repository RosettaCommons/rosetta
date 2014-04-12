// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/full_model.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_full_model_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_full_model_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace full_model { extern BooleanOptionKey const full_model; }
namespace full_model { extern IntegerVectorOptionKey const cutpoint_open; }
namespace full_model { extern IntegerVectorOptionKey const cutpoint_closed; }
namespace full_model { extern StringVectorOptionKey const other_poses; }
namespace full_model { extern IntegerVectorOptionKey const extra_min_res; }
namespace full_model { extern IntegerVectorOptionKey const jump_res; }
namespace full_model { extern IntegerVectorOptionKey const root_res; }
namespace full_model { extern IntegerVectorOptionKey const virtual_sugar_res; }
namespace full_model { extern IntegerVectorOptionKey const virtual_res; }
namespace full_model { extern IntegerVectorOptionKey const sample_res; }
namespace full_model { extern IntegerVectorOptionKey const calc_rms_res; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
