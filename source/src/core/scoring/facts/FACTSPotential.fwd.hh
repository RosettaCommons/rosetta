// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @file:   core/scoring/facts/FACTSPotential.fwd.hh
// @brief:  This file only creats some easy to use aliases
// @author: Hahnbeom Park


#ifndef INCLUDED_core_scoring_facts_FACTSPotential_FWD_HH
#define INCLUDED_core_scoring_facts_FACTSPotential_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <core/chemical/ResidueType.hh>

namespace core {
namespace scoring {

class FACTSRsdTypeMap;

//Declaring a class of type FACTSResidueInfo
class FACTSRsdTypeInfo;
//Creating an alias for a pointer of type FACTSResidueInfo
typedef utility::pointer::shared_ptr< FACTSRsdTypeInfo > FACTSRsdTypeInfoOP;
typedef utility::pointer::shared_ptr< FACTSRsdTypeInfo const > FACTSRsdTypeInfoCOP;

//Declaring a class of type FACTSResidueInfo
class FACTSResidueInfo;
//Creating an alias for a pointer of type FACTSResidueInfo
typedef utility::pointer::shared_ptr< FACTSResidueInfo > FACTSResidueInfoOP;

//Declaring a class of type FACTSPoseInfo
class FACTSPoseInfo;
//Creating an alias for a pointer of type FACTSPoseInfo
typedef utility::pointer::shared_ptr< FACTSPoseInfo > FACTSPoseInfoOP;

class FACTSRotamerSetInfo;
typedef utility::pointer::shared_ptr< FACTSRotamerSetInfo > FACTSRotamerSetInfoOP;

//Declaring a class of type FACTSPotential
class FACTSPotential;
//Creating an alias for a pointer of type FACTSResidueInfo
typedef utility::pointer::shared_ptr< FACTSPotential > FACTSPotentialOP;

} // scoring
} // core

#endif
