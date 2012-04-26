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
#include <core/pack/task/operation/TaskOperationRegistrator.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/moves/MoverFactory.hh>

//mover creators
#include <devel/matdes/SymmetrizerMoverCreator.hh>
#include <devel/matdes/TaskAwareSymMinMoverCreator.hh>

// Utility Headers

// Task Operation creators
#include <devel/matdes/BuildingBlockInterfaceOperationCreator.hh>
#include <devel/matdes/RestrictToNonzeroSASAOperationCreator.hh>
#include <devel/matdes/RestrictIdentitiesToRepackingOperationCreator.hh>

// Filter creators
#include <devel/matdes/OligomericAverageDegreeFilterCreator.hh>

#include <utility/vector1.hh>


namespace devel {

// Mover creators
static protocols::moves::MoverRegistrator< devel::matdes::SymmetrizerMoverCreator > reg_SymmetrizerMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::TaskAwareSymMinMoverCreator > reg_TaskAwareSymMinMoverCreator;

// Task creators
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::BuildingBlockInterfaceOperationCreator > BuildingBlockInterfaceOperationCreator_registrator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::RestrictToNonzeroSASAOperationCreator > RestrictToNonzeroSASAOperationCreator_registrator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::RestrictIdentitiesToRepackingOperationCreator > RestrictIdentitiesToRepackingOperationCreator_registrator;

// Filter creators
static protocols::filters::FilterRegistrator< devel::matdes::OligomericAverageDegreeFilterCreator > OligomericAverageDegreeFilterCreator_registrator;

void init( int argc, char * argv [] )
{
	protocols::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	protocols::init::init( args );
} // init

} // devel

