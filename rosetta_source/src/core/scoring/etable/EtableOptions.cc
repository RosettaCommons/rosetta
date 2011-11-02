// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <core/scoring/etable/EtableOptions.hh>

// AUTO-REMOVED #include <core/scoring/types.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace etable {

EtableOptions::EtableOptions() :
	max_dis( 6.0 ),
	bins_per_A2( 20 ),
	Wradius( 1.0 ),
	lj_switch_dis2sigma( 0.6 ),
	disable_polar_desolvation( false )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	max_dis = option[ score::fa_max_dis ];
	if( option[ score::no_smooth_etables ] && !option[ score::fa_max_dis ].user() ) {
		basic::T("core.scoring.etable") << "no_smooth_etables requested and fa_max_dis not specified: using 5.5 as default" << std::endl;
		max_dis = 5.5;
	}
	if ( option[ score::no_lk_polar_desolvation ] ) {
		disable_polar_desolvation = false;
	}
}

} // etable
} // scoring
} // core


