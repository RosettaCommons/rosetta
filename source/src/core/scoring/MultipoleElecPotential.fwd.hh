// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/MultipoleElecPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_MultipoleElecPotential_fwd_hh
#define INCLUDED_core_scoring_MultipoleElecPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

///
class MultipoleElecResidueInfo;
typedef utility::pointer::shared_ptr< MultipoleElecResidueInfo > MultipoleElecResidueInfoOP;
typedef utility::pointer::shared_ptr< MultipoleElecResidueInfo const > MultipoleElecResidueInfoCOP;

///
class MultipoleElecPoseInfo;
typedef utility::pointer::shared_ptr< MultipoleElecPoseInfo > MultipoleElecPoseInfoOP;
typedef utility::pointer::shared_ptr< MultipoleElecPoseInfo const > MultipoleElecPoseInfoCOP;

///
class MultipoleElecRotamerSetInfo;
typedef utility::pointer::shared_ptr< MultipoleElecRotamerSetInfo > MultipoleElecRotamerSetInfoOP;

///
class MultipoleParameter;
typedef utility::pointer::shared_ptr< MultipoleParameter > MultipoleParameterOP;

///
class MultipoleElecPotential;

typedef  utility::pointer::shared_ptr< MultipoleElecPotential > MultipoleElecPotentialOP;
typedef  utility::pointer::shared_ptr< MultipoleElecPotential const > MultipoleElecPotentialCOP;



} // scoring
} // core

#endif
