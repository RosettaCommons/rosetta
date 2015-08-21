// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_evaluation_PCA_fwd_hh
#define INCLUDED_protocols_evaluation_PCA_fwd_hh


// Unit Headers
// #include <protocols/evaluation/PCA.fwd.hh>

// Package Headers

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace protocols {
namespace evaluation {

class PCA;
typedef utility::pointer::shared_ptr< PCA > PCA_OP;

}
}

#endif
