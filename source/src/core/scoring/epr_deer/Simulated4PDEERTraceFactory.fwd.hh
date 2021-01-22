// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/Simulated4PDEERTraceFactory.cc
/// @brief  Class that simulates the observable signal resulting from electron dipolar coupling
/// @details The dipolar coupling signal between two electrons can be simulated from a distance
///       distribution. This class does that simulation by pulling values from the datacache
///       container objects (DEERData) and effectively storing a vector with values of interest.
/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_fwd_hh
#define INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace epr_deer {

class Simulated4PDEERTraceFactory;

typedef utility::pointer::shared_ptr< Simulated4PDEERTraceFactory > Simulated4PDEERTraceFactoryOP;
typedef utility::pointer::shared_ptr< Simulated4PDEERTraceFactory const > Simulated4PDEERTraceFactoryCOP;


} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_epr_deer_Simulated4PDEERTraceFactory_fwd_hh
