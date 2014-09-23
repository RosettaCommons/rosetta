// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MMLJEnergyInter.fwd.hh
/// @brief  molecular mechanics lj energy forward declaration
/// @author P. Douglas Renfrew (renfrew@nyu.edu)

#ifndef INCLUDED_core_scoring_methods_MMLJEnergyInter_fwd_hh
#define INCLUDED_core_scoring_methods_MMLJEnergyInter_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class MMLJEnergyInter;

typedef  utility::pointer::weak_ptr< MMLJEnergyInter > MMLJEnergyInterAP;
typedef  utility::pointer::weak_ptr< MMLJEnergyInter const > MMLJEnergyInterCAP;
typedef  utility::pointer::shared_ptr< MMLJEnergyInter > MMLJEnergyInterOP;
typedef  utility::pointer::shared_ptr< MMLJEnergyInter const > MMLJEnergyInterCOP;

} // namespace methods
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_methods_MMLJEnergyInter_FWD_HH
