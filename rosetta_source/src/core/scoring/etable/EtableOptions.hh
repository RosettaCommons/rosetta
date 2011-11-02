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

#ifndef INCLUDED_core_scoring_etable_EtableOptions_hh
#define INCLUDED_core_scoring_etable_EtableOptions_hh

// AUTO-REMOVED #include <core/scoring/types.hh>

#include <core/types.hh>


namespace core {
namespace scoring {
namespace etable {

class EtableOptions {

public:

	EtableOptions();

public:

	Real max_dis;
	int bins_per_A2;
	Real Wradius;
	Real lj_switch_dis2sigma;
	bool disable_polar_desolvation;

private:

};

} // etable
} // scoring
} // core

#endif
