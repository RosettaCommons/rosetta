// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/VdWTinkerPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_VdWTinkerPotential_fwd_hh
#define INCLUDED_core_scoring_VdWTinkerPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

///
class VdWTinkerResidueInfo;
typedef utility::pointer::shared_ptr< VdWTinkerResidueInfo > VdWTinkerResidueInfoOP;
typedef utility::pointer::shared_ptr< VdWTinkerResidueInfo const > VdWTinkerResidueInfoCOP;

///
class VdWTinkerPoseInfo;
typedef utility::pointer::shared_ptr< VdWTinkerPoseInfo > VdWTinkerPoseInfoOP;
typedef utility::pointer::shared_ptr< VdWTinkerPoseInfo const > VdWTinkerPoseInfoCOP;

///
class VdWTinkerRotamerSetInfo;
typedef utility::pointer::shared_ptr< VdWTinkerRotamerSetInfo > VdWTinkerRotamerSetInfoOP;

///
class MultipoleParameter;
typedef utility::pointer::shared_ptr< MultipoleParameter > MultipoleParameterOP;

///
class VdWTinkerPotential;

typedef  utility::pointer::shared_ptr< VdWTinkerPotential > VdWTinkerPotentialOP;
typedef  utility::pointer::shared_ptr< VdWTinkerPotential const > VdWTinkerPotentialCOP;



} // scoring
} // core

#endif
