// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/custom_pair_distance/FullatomCustomPairDistanceEnergy.fwd.hh
/// @brief  class declaration
/// @author David E Kim


#ifndef INCLUDED_core_scoring_custom_pair_distance_FullatomCustomPairDistanceEnergy_fwd_hh
#define INCLUDED_core_scoring_custom_pair_distance_FullatomCustomPairDistanceEnergy_fwd_hh

//Utility Headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace custom_pair_distance {

struct custom_pair_and_func_struct;
struct resatom_and_func_struct;

class FullatomCustomPairDistanceEnergy;

typedef utility::pointer::shared_ptr< FullatomCustomPairDistanceEnergy > FullatomCustomPairDistanceEnergyOP;
typedef utility::pointer::shared_ptr< FullatomCustomPairDistanceEnergy const > FullatomCustomPairDistanceEnergyCOP;

class DistanceFunc;

typedef utility::pointer::shared_ptr< DistanceFunc > DistanceFuncOP;
typedef utility::pointer::shared_ptr< DistanceFunc const > DistanceFuncCOP;
typedef utility::pointer::weak_ptr< DistanceFunc > DistanceFuncAP;
typedef utility::pointer::weak_ptr< DistanceFunc const > DistanceFuncCAP;


class CacheableAtomPairFuncMap;

typedef utility::pointer::shared_ptr< CacheableAtomPairFuncMap > CacheableAtomPairFuncMapOP;
typedef utility::pointer::shared_ptr< CacheableAtomPairFuncMap const > CacheableAtomPairFuncMapCOP;

}
}
}

#endif
