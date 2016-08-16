// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/GenBornPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_GenBornPotential_fwd_hh
#define INCLUDED_core_scoring_GenBornPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {


class GenBornResidueInfo;
typedef utility::pointer::shared_ptr< GenBornResidueInfo > GenBornResidueInfoOP;


class GenBornPoseInfo;
typedef utility::pointer::shared_ptr< GenBornPoseInfo > GenBornPoseInfoOP;


class GenBornRotamerSetInfo;
typedef utility::pointer::shared_ptr< GenBornRotamerSetInfo > GenBornRotamerSetInfoOP;


class GenBornPotential;

typedef  utility::pointer::shared_ptr< GenBornPotential > GenBornPotentialOP;
typedef  utility::pointer::shared_ptr< GenBornPotential const > GenBornPotentialCOP;


} // scoring
} // core

#endif
