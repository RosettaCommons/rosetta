// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/StubStubType.fwd.hh
/// @brief For computing or constrainign rotation/translation -- different stub-stube conventions in use.
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_StubStubType_FWD_HH
#define INCLUDED_core_pose_rna_StubStubType_FWD_HH

namespace core {
namespace pose {
namespace rna {

/// @details
// BASE_CENTROID   Rosetta standard base coordinate systems.
/// O3P_TO_O5P      ...-C4'-C3'-O3' to O5'-C5'-C4'-... Centers are O3' and O5'.
/// CHAINBREAK      O3'-OVL1-OVL2 to OVU-P-O5', centered on OVL1 and P [overlap means transl = 0, rotation = I ]
///
enum StubStubType {
	NONE,
	BASE_CENTROID,
	O3P_TO_O5P,
	CHAINBREAK
};

} //rna
} //pose
} //core

#endif
