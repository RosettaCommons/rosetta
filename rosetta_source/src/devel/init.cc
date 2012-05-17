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
#include <protocols/jd2/parser/DataLoaderFactory.hh>
//#include <devel/constrained_sequence_design/SequenceConstraintFactory.hh>

//mover creators
#include <devel/enzdes/EnzdesRemodelMoverCreator.hh>
//#include <devel/constrained_sequence_design/ConstrainedDesignMoverCreator.hh>
#include <devel/matdes/SymmetrizerMoverCreator.hh>
#include <devel/matdes/TaskAwareSymMinMoverCreator.hh>
#include <devel/matdes/StoreTaskMoverCreator.hh>

// Filter creators
#include <devel/matdes/OligomericAverageDegreeFilterCreator.hh>
#include <devel/matdes/SymUnsatHbondFilterCreator.hh>
#include <devel/matdes/AverageInterfaceEnergyFilterCreator.hh>
#include <devel/matdes/TaskAwareAlaScanCreator.hh>
#include <devel/matdes/SaveResfileToDiskFilterCreator.hh>
#include <devel/matdes/TaskAwareSASAFilterCreator.hh>
#include <devel/matdes/InterfacePackingFilterCreator.hh>

// dataloader creators
//#include <devel/constrained_sequence_design/SequenceConstraintLoaderCreator.hh>

// SequenceConstraint creators
//#include <devel/constrained_sequence_design/constraints/MaximunNumberPerResidueTypeConstraintCreator.hh>

// Utility Headers

// Task Operation creators
#include <devel/znhash/SymmZnMoversAndTaskOpsCreators.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>
#include <devel/matdes/BuildingBlockInterfaceOperationCreator.hh>
#include <devel/matdes/RestrictToNonzeroSASAOperationCreator.hh>
#include <devel/matdes/RestrictIdentitiesOperationCreator.hh>
#include <devel/matdes/RetrieveStoredTaskOperationCreator.hh>

#include <utility/vector1.hh>


namespace devel {

// Mover creators
protocols::moves::MoverRegistrator< enzdes::EnzdesRemodelMoverCreator > reg_EnzdesRemodelMoverCreator;
protocols::moves::MoverRegistrator< vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator > reg_LoadVarSolDistSasaCalculatorMoverCreator;
protocols::moves::MoverRegistrator< devel::znhash::InsertZincCoordinationRemarkLinesCreator > reg_InsertZincCoordinationRemarkLinesCreator;
protocols::moves::MoverRegistrator< znhash::LoadZnCoordNumHbondCalculatorMoverCreator > reg_LoadZnCoordNumHbondCalculatorMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::SymmetrizerMoverCreator > reg_SymmetrizerMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::TaskAwareSymMinMoverCreator > reg_TaskAwareSymMinMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::StoreTaskMoverCreator > reg_StoreTaskMoverCreator;

// Task creators
core::pack::task::operation::TaskOperationRegistrator< devel::znhash::DisableZnCoordinationResiduesTaskOpCreator > reg_DisableZnCoordinationResiduesTaskOpCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::BuildingBlockInterfaceOperationCreator > BuildingBlockInterfaceOperationCreator_registrator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::RestrictToNonzeroSASAOperationCreator > RestrictToNonzeroSASAOperationCreator_registrator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::RestrictIdentitiesOperationCreator > RestrictIdentitiesOperationCreator_registrator;
static core::pack::task::operation::TaskOperationRegistrator< devel::matdes::RetrieveStoredTaskOperationCreator > RetrieveStoredTaskOperationCreator_registrator;
 
// Filter creators 
static protocols::filters::FilterRegistrator< devel::matdes::OligomericAverageDegreeFilterCreator > OligomericAverageDegreeFilterCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::SymUnsatHbondFilterCreator > SymUnsatHbondFilterCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::AverageInterfaceEnergyFilterCreator > AverageInterfaceEnergyFilterCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::TaskAwareAlaScanCreator > TaskAwareAlaScanCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::SaveResfileToDiskFilterCreator > SaveResfileToDiskFilterCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::TaskAwareSASAFilterCreator > TaskAwareSASAFilterCreator_registrator;
static protocols::filters::FilterRegistrator< devel::matdes::InterfacePackingFilterCreator > InterfacePackingFilterCreator_registrator;

void init( int argc, char * argv [] )
{
	protocols::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	protocols::init::init( args );
} // init

} // devel

