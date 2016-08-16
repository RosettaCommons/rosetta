// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/etable/EtableOptions.fwd.hh
/// @brief forward header for EtableOptions class
//  @author

#ifndef INCLUDED_core_scoring_etable_EtableOptions_fwd_hh
#define INCLUDED_core_scoring_etable_EtableOptions_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace etable {

class EtableOptions;

typedef utility::pointer::shared_ptr< EtableOptions > EtableOptionsOP;
typedef utility::pointer::shared_ptr< EtableOptions const > EtableOptionsCOP;

} //etable
} //scoring
} //core

#endif // INCLUDED_core_scoring_etable_EtableOptions_FWD_HH
