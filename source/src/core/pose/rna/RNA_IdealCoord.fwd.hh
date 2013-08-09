// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rna/RNA_IdealCoord.fwd.hh
/// @brief Apply ideal RNA geometry to a residue or a pose
/// @author Fang-Chieh Chou


#ifndef INCLUDED_core_pose_rna_RNA_IdealCoord_fwd_HH
#define INCLUDED_core_pose_rna_RNA_IdealCoord_fwd_HH

// C++

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.fwd.hh>

namespace core {
namespace pose {
namespace rna {

class RNA_IdealCoord;
typedef utility::pointer::owning_ptr< RNA_IdealCoord > RNA_IdealCoordOP;
typedef utility::pointer::access_ptr< RNA_IdealCoord > RNA_IdealCoordAP;

}
}
}

#endif
