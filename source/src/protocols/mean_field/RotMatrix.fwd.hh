// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RotMatrix.fwd.hh

/// @brief   Forward declarations for RotMatrix
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_RotMatrix_fwd_hh
#define INCLUDED_protocols_mean_field_RotMatrix_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace mean_field {

/// @brief Forward declaration for RotMatrix, a jagged array of RotProbs used for mean-field calculations
class RotMatrix;
typedef utility::pointer::shared_ptr< RotMatrix > RotMatrixOP;
typedef utility::pointer::shared_ptr< RotMatrix const > RotMatrixCOP;

}  //namespace mean_field
}  //namespace protocols

#endif  // INCLUDED_protocols_mean_field_RotMatrix_fwd_hh
