// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/SheetConstraintsRCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009

// Unit header
#include <protocols/forge/constraints/SheetConstraintsRCG.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/fldsgn/BluePrint.hh>
#include <protocols/flxbb/utility.hh>
#include <core/id/SequenceMapping.hh>

// Project headers
#include <basic/Tracer.hh>

#include <protocols/fldsgn/topology/HSSTriplet.hh>

static basic::Tracer TR( "protocols.forge.constraints.SheetConstraintsRCG" );

namespace protocols{
namespace forge{
namespace constraints{

/// @brief
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue ):
	RemodelConstraintGenerator(),
	blueprint_( blue ),
	coef_( 1.0 ),
	dist_( 5.5 )
{}

/// @brief value constructor
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue, Real const coef ):
	RemodelConstraintGenerator(),
	blueprint_( blue ),
	coef_( coef ),
	dist_( 5.5 )
{}

/// @brief value constructor
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue, Real const coef, Real const dist ):
	RemodelConstraintGenerator(),
	blueprint_( blue ),
	coef_( coef ),
	dist_( dist )
{}

/// @brief
SheetConstraintsRCG::~SheetConstraintsRCG() {}

/// @brief
void
SheetConstraintsRCG::set_blueprint( BluePrintOP const & blue )
{
	blueprint_ = blue;
}

/// @brief set weight
void
SheetConstraintsRCG::set_weight( Real const coef )
{
	coef_ = coef;
}

/// @brief set distance of constraint
void
SheetConstraintsRCG::set_distance( Real const dist )
{
	dist_ = dist;
}

/// @brief
void
SheetConstraintsRCG::generate_remodel_constraints( Pose const & pose )
{
	using core::scoring::constraints::ConstraintOPs;
	ConstraintOPs constraints = protocols::flxbb::constraints_sheet( pose, blueprint_, coef_, dist_ );
	this->add_constraints( constraints );
} //generate constraints


} //namespace constraints
} //namespace forge
} //namespace protocols
