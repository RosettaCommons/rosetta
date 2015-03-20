// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/util.hh
/// @brief util functions for building structures from components
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_util_hh
#define INCLUDED_protocols_denovo_design_util_hh

// Unit headers

// Protocol headers

// Package headers

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Numeric/Utility Headers

// C++ Headers

namespace protocols {
namespace denovo_design {

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
void construct_poly_ala_pose( core::pose::Pose & pose, bool const keep_disulf );

} // denovo_design
} // protocols

#endif
