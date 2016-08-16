// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/compute_sasa.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_packing_score_params_hh
#define INCLUDED_core_scoring_packstat_packing_score_params_hh

#include <core/scoring/packstat/PackingScore.hh>

namespace core {
namespace scoring {
namespace packstat {

void init_packing_score_respred( PackingScore & p );
void init_packing_score_discrim( PackingScore & p );

} // namesapce packstat
} // namespace scoring
} // namespace core

#endif
