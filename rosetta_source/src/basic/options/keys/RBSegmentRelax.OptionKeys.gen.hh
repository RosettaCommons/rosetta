// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_RBSegmentRelax_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_RBSegmentRelax_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace RBSegmentRelax { extern BooleanOptionKey const RBSegmentRelax; }
namespace RBSegmentRelax { extern FileOptionKey const input_pdb; }
namespace RBSegmentRelax { extern FileOptionKey const rb_file; }
namespace RBSegmentRelax { extern RealOptionKey const cst_wt; }
namespace RBSegmentRelax { extern RealOptionKey const cst_width; }
namespace RBSegmentRelax { extern StringOptionKey const cst_pdb; }
namespace RBSegmentRelax { extern IntegerOptionKey const nrbmoves; }
namespace RBSegmentRelax { extern IntegerOptionKey const nrboutercycles; }
namespace RBSegmentRelax { extern StringOptionKey const rb_scorefxn; }
namespace RBSegmentRelax { extern BooleanOptionKey const skip_fragment_moves; }
namespace RBSegmentRelax { extern BooleanOptionKey const skip_seqshift_moves; }
namespace RBSegmentRelax { extern BooleanOptionKey const skip_rb_moves; }
namespace RBSegmentRelax { extern RealVectorOptionKey const helical_movement_params; }
namespace RBSegmentRelax { extern RealVectorOptionKey const strand_movement_params; }
namespace RBSegmentRelax { extern RealVectorOptionKey const default_movement_params; }
namespace RBSegmentRelax { extern IntegerOptionKey const cst_seqwidth; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
