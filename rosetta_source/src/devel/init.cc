// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/init.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <devel/init.hh>
#include <protocols/init/init.hh>

// Factories
#include <protocols/moves/MoverFactory.hh>
#include <protocols/jd2/parser/DataLoaderFactory.hh>
//#include <devel/constrained_sequence_design/SequenceConstraintFactory.hh>

//mover creators
#include <devel/enzdes/EnzdesRemodelMoverCreator.hh>
//#include <devel/constrained_sequence_design/ConstrainedDesignMoverCreator.hh>

// dataloader creators
//#include <devel/constrained_sequence_design/SequenceConstraintLoaderCreator.hh>

// SequenceConstraint creators
//#include <devel/constrained_sequence_design/constraints/MaximunNumberPerResidueTypeConstraintCreator.hh>

// Utility Headers
#include <utility/vector1.hh>


namespace devel {

// movers
static protocols::moves::MoverRegistrator< enzdes::EnzdesRemodelMoverCreator > reg_EnzdesRemodelMoverCreator;
//static protocols::moves::MoverRegistrator< constrained_sequence_design::ConstrainedDesignMoverCreator > reg_ConstrainedDesignMoverCreator;

// SequenceConstraints
//static constrained_sequence_design::SequenceConstraintRegistrator< constrained_sequence_design::constraints::MaximunNumberPerResidueTypeConstraintCreator > reg_MaximunNumberPerResidueTypeConstraint;

// data loaders
//static protocols::jd2::parser::DataLoaderRegistrator< constrained_sequence_design::SequenceConstraintLoaderCreator > reg_SequenceConstraintLoaderCreator;

void init( int argc, char * argv [] )
{
	protocols::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	protocols::init::init( args );
} // init

} // devel

