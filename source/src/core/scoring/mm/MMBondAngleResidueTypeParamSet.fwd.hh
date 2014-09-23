// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh
/// @brief  core::scoring::mm::MMBondAngleResidueTypeParamSet forward declarations
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParamSet_fwd_hh
#define INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParamSet_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace scoring {
namespace mm {


// Forward
class MMBondAngleResidueTypeParamSet;

typedef  utility::pointer::weak_ptr< MMBondAngleResidueTypeParamSet > MMBondAngleResidueTypeParamSetAP;
typedef  utility::pointer::weak_ptr< MMBondAngleResidueTypeParamSet const > MMBondAngleResidueTypeParamSetCAP;
typedef  utility::pointer::shared_ptr< MMBondAngleResidueTypeParamSet > MMBondAngleResidueTypeParamSetOP;
typedef  utility::pointer::shared_ptr< MMBondAngleResidueTypeParamSet const > MMBondAngleResidueTypeParamSetCOP;

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_MMBondAngleResidueTypeParamSet_FWD_HH
