// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/chemical/carbohydrates/CarbohydrateInfo.fwd.hh
/// @brief   Forward declarations for CarbohydrateInfo.
/// @author  Labonte

#ifndef INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_FWD_HH
#define INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace carbohydrates {

/// @brief  A class containing carbohydrate-specific information.
class CarbohydrateInfo;

typedef utility::pointer::shared_ptr< CarbohydrateInfo > CarbohydrateInfoOP;
typedef utility::pointer::shared_ptr< CarbohydrateInfo const > CarbohydrateInfoCOP;

}  // namespace carbohydrates
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_carbohydrates_CarbohydrateInfo_FWD_HH
