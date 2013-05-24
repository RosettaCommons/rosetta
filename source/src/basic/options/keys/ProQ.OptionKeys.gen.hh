// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/ProQ.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_ProQ_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_ProQ_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace ProQ { extern BooleanOptionKey const ProQ; }
namespace ProQ { extern IntegerOptionKey const svmmodel; }
namespace ProQ { extern StringOptionKey const basename; }
namespace ProQ { extern BooleanOptionKey const membrane; }
namespace ProQ { extern BooleanOptionKey const prof_bug; }
namespace ProQ { extern BooleanOptionKey const output_feature_vector; }
namespace ProQ { extern BooleanOptionKey const output_local_prediction; }
namespace ProQ { extern StringOptionKey const prefix; }
namespace ProQ { extern BooleanOptionKey const use_gzip; }
namespace ProQ { extern RealOptionKey const normalize; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
