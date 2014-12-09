// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_SRFD_Implementations_HH
#define INCLUDED_core_fragment_SRFD_Implementations_HH

// Unit Headers
#include <core/fragment/SingleResidueFragData.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh> // for ResidueSRFD

#include <core/kinematics/types.hh>
#include <core/id/TorsionID.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace fragment {

/// FURTHER EXAMPLES FOR SRFDs ... ///
/*
// not used right now
class ResidueSRFD : public SingleResidueFragData {
public:
	bool apply( pose::Pose&, Size seq_pos );
private:
	conformation::Residue data_;
};

// very specific
class BBProteinTorsionSRFD : public SingleResidueFragData {
public:

private:
	Real phi_;
	Real psi_;
	Real omega_;
	// should it contain ss-type ? how do we want to handle this now?
	char secstruct_;
};


class RNATorsionSRFD : public BaseTorsionSRFD {
	//???
};
*/

///class GeneralDofSRFD : public SingleResidueFragData {
	///* similar to the TorsionSRFD but use also the DOF_TYPE */
	///* CAREFUL: this has to be implemented such that there is no dependence on atom-tree
	//	 */
	/* this, however, might need expansion of the pose interface
			 enum DOF_Type {
			 PHI = 1, // used for lookup into utility::vector1
			 THETA,
			 D,
			 RB1,
			 RB2,
			 RB3,
			 RB4,
			 RB5,
			 RB6
			 }; */
//};

} //fragment
} //core

#endif
