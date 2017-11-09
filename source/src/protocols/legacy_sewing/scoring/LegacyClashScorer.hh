// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyClashScorer.hh
///
/// @brief Identifies and scores backbone clashes in LEGACY_SEWING Assemblies
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_scoring_LegacyClashScorer_hh
#define INCLUDED_protocols_legacy_sewing_scoring_LegacyClashScorer_hh

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyClashScorer.fwd.hh>
#include <protocols/legacy_sewing/scoring/LegacyAssemblyScorer.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Utility headers
#include <utility/vector1.hh>

#include <utility/tag/Tag.hh>

namespace protocols {
namespace legacy_sewing  {
namespace scoring {

class LegacyClashScorer : public LegacyAssemblyScorer {

public:

	///@brief default construct
	LegacyClashScorer();

	virtual ~LegacyClashScorer(){}

	virtual
	core::Real
	score(
		AssemblyCOP assembly
	);
};


} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
