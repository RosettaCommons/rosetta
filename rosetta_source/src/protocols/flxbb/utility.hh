// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/utility.hh
/// @brief
/// @author Nobuyasu Koga (nobuyasu@uw.edu)

#ifndef INCLUDED_protocols_flxbb_utility_hh
#define INCLUDED_protocols_flxbb_utility_hh

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

#include <protocols/fldsgn/BluePrint.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flxbb {

typedef core::Real Real;
typedef core::pose::Pose Pose;
typedef protocols::fldsgn::BluePrintOP BluePrintOP;
typedef core::scoring::constraints::ConstraintOPs ConstraintOPs;

/// @brief constrain hydrogen bonds in beta sheet  ( not used )
//void constraints_sheet( pose::Pose & pose, scoring::ScoreFunctionOP & scorefxn, BluePrintOP & blueprint_ );

/// @brief constrain between Ca atoms in beta sheet, which are specified in blueprint file
ConstraintOPs
constraints_sheet( Pose const & pose, BluePrintOP const & blue, Real const coef, Real const condist=5.5 );

/// @brief constrain between Ca atoms in beta sheet
ConstraintOPs
constraints_sheet( Pose const & pose, Real const coef, Real const condist=5.5 );

/// @brief constraint between N- and C-terminal Ca atoms
ConstraintOPs
constraints_NtoC( Pose const & pose, Real const coef, Real const condist=11.0 );


} // namespace flxbb
} // namespace protocols

#endif
