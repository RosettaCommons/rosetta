// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/boinc.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_boinc_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_boinc_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace boinc { extern BooleanOptionKey const boinc; }
namespace boinc { extern BooleanOptionKey const graphics; }
namespace boinc { extern BooleanOptionKey const fullscreen; }
namespace boinc { extern IntegerOptionKey const max_fps; }
namespace boinc { extern IntegerOptionKey const max_cpu; }
namespace boinc { extern BooleanOptionKey const noshmem; }
namespace boinc { extern IntegerOptionKey const cpu_run_time; }
namespace boinc { extern IntegerOptionKey const max_nstruct; }
namespace boinc { extern RealOptionKey const cpu_frac; }
namespace boinc { extern RealOptionKey const frame_rate; }
namespace boinc { extern BooleanOptionKey const watchdog; }
namespace boinc { extern IntegerOptionKey const watchdog_time; }
namespace boinc { extern IntegerOptionKey const cpu_run_timeout; }
namespace boinc { extern FileOptionKey const description_file; }
namespace boinc { extern RealOptionKey const score_cut_pct; }
namespace boinc { extern FileOptionKey const score_cut_fl; }
namespace boinc { extern BooleanOptionKey const score_cut_smart_throttle; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
