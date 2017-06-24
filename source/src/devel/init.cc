// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/devel/init.cc
/// @brief  Declare WidgetRegistrators as static (global) variables in this .cc file
///         so that at load time, they will be initialized, and the Creator classes
///         they register will be handed to the appropriate WidgetFactory.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <devel/init.hh>
#include <devel/svn_version.hh>
#include <protocols/init/init.hh>

// Factories
#include <core/pack/task/operation/TaskOperationRegistrator.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/scoring/methods/EnergyMethodRegistrator.hh>
#include <protocols/filters/FilterFactory.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/parser/DataLoaderFactory.hh>
//#include <devel/constrained_sequence_design/SequenceConstraintFactory.hh>

//mover creators
#include <devel/denovo_design/ConnectJumpsCreator.hh>
#include <devel/denovo_design/DumpStatsSSCreator.hh>
#include <devel/denovo_design/RestrictRegionCreator.hh>
#include <devel/domain_insertion/FusePosesNtoCMoverCreator.hh>
#include <devel/domain_insertion/InsertionSiteTestMoverCreator.hh>
#include <devel/enzdes/EnzdesRemodelMoverCreator.hh>
#include <devel/loop_creation/LoopCreationMoverCreator.hh>
#include <devel/loop_creation/LoophashLoopInserterCreator.hh>
#include <devel/loop_creation/IterativeLoophashLoopInserterCreator.hh>
#include <devel/loop_creation/CCDLoopCloserCreator.hh>
#include <devel/loop_creation/FragmentLoopInserterCreator.hh>
#include <devel/matdes/SymmetrizerMoverCreator.hh>
#include <devel/replica_docking/AddEncounterConstraintMoverCreator.hh>
#include <devel/replica_docking/ModulatedMoverCreator.hh>
#include <devel/loophash_loopclosure/LoopHashLoopClosureMoverCreator.hh>
#include <devel/splice/SpliceCreator.hh>
#include <devel/splice/DesignInterfacesOperationCreator.hh>
#include <devel/splice/FindEndpointsOperationCreator.hh>
#include <devel/splice/AlignEndsMoverCreator.hh>
#include <devel/cutoutdomain/CutOutDomainCreator.hh>
#include <devel/splice/RBOutMoverCreator.hh>
#include <devel/splice/RBInMoverCreator.hh>
#include <devel/matdes/StoreQuasiSymmetricTaskMoverCreator.hh>

// Filter creators
#include <devel/denovo_design/filters/CavityVolumeFilterCreator.hh>
#include <devel/denovo_design/filters/CoreResiduesPerElementFilterCreator.hh>
#include <devel/denovo_design/filters/FoldabilityFilterCreator.hh>
#include <devel/denovo_design/filters/SSShapeComplementarityFilterCreator.hh>
#include <devel/replica_docking/InteractionScoreFilterCreator.hh>
#include <devel/matdes/GenericSymmetricSamplerCreator.hh>
#include <devel/replica_docking/IrmsdFilterCreator.hh>
#include <devel/replica_docking/CaIrmsdFilterCreator.hh>
#include <devel/buns/BuriedUnsatHbondFilter2Creator.hh>
#include <devel/replica_docking/FnatFilterCreator.hh>
#include <devel/replica_docking/LrmsdFilterCreator.hh>
#include <devel/replica_docking/FnonnatFilterCreator.hh>
#include <devel/replica_docking/WrapFilterAsEvaluatorCreator.hh>

// Energy method creators
#include <devel/denovo_design/scoring/SideChainNeighborsEnergyCreator.hh>

// dataloader creators
//#include <devel/constrained_sequence_design/SequenceConstraintLoaderCreator.hh>

// SequenceConstraint creators
//#include <devel/constrained_sequence_design/constraints/MaximunNumberPerResidueTypeConstraintCreator.hh>

// Utility Headers

// Task Operation creators
#include <devel/denovo_design/task_operations/HighestEnergyRegionCreator.hh>
#include <devel/znhash/SymmZnMoversAndTaskOpsCreators.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>

#include <utility/vector1.hh>

namespace devel {

// Mover creators
static protocols::moves::MoverRegistrator< denovo_design::ConnectJumpsCreator > reg_ConnectJumpsCreator;
static protocols::moves::MoverRegistrator< denovo_design::DumpStatsSSCreator > reg_DumpStatsSSCreator;
static protocols::moves::MoverRegistrator< denovo_design::RestrictRegionCreator > reg_RestrictRegionCreator;
protocols::moves::MoverRegistrator< domain_insertion::FusePosesNtoCMoverCreator > reg_FusePosesNtoCMoverCreator;
protocols::moves::MoverRegistrator< domain_insertion::InsertionSiteTestMoverCreator > reg_InsertionSiteTestMoverCreator;
protocols::moves::MoverRegistrator< domain_insertion::SetupCoiledCoilFoldTreeMoverCreator > reg_SetupCoiledCoilFoldTreeMoverCreator;
protocols::moves::MoverRegistrator< enzdes::EnzdesRemodelMoverCreator > reg_EnzdesRemodelMoverCreator;
protocols::moves::MoverRegistrator< vardist_solaccess::LoadVarSolDistSasaCalculatorMoverCreator > reg_LoadVarSolDistSasaCalculatorMoverCreator;
protocols::moves::MoverRegistrator< devel::znhash::InsertZincCoordinationRemarkLinesCreator > reg_InsertZincCoordinationRemarkLinesCreator;
protocols::moves::MoverRegistrator< znhash::LoadZnCoordNumHbondCalculatorMoverCreator > reg_LoadZnCoordNumHbondCalculatorMoverCreator;
static protocols::moves::MoverRegistrator< devel::loop_creation::LoopCreationMoverCreator > reg_LoopCreationMoverCreator;
static protocols::moves::MoverRegistrator< devel::loop_creation::FragmentLoopInserterCreator > reg_FragmentLoopInserterCreator;
static protocols::moves::MoverRegistrator< devel::loop_creation::CCDLoopCloserCreator > reg_CCDLoopCloserCreator;
static protocols::moves::MoverRegistrator< devel::loop_creation::LoophashLoopInserterCreator > reg_LoophashLoopInserterCreator;
static protocols::moves::MoverRegistrator< devel::loop_creation::IterativeLoophashLoopInserterCreator > reg_IterativeLoophashLoopInserterCreator;
static protocols::moves::MoverRegistrator< devel::matdes::SymmetrizerMoverCreator > reg_SymmetrizerMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::StoreQuasiSymmetricTaskMoverCreator > reg_StoreQuasiSymmetricTaskMoverCreator;
static protocols::moves::MoverRegistrator< devel::matdes::GenericSymmetricSamplerCreator > reg_GenericSymmetricSamplerCreator;
static protocols::moves::MoverRegistrator< replica_docking::AddEncounterConstraintMoverCreator > reg_AddEncounterConstraintMoverCreator;
static protocols::moves::MoverRegistrator< replica_docking::ModulatedMoverCreator > reg_ModulatedMoverCreator;
static protocols::moves::MoverRegistrator< loophash_loopclosure::LoopHashLoopClosureMoverCreator > reg_LoopHashLoopClosureMoverCreator;
static protocols::moves::MoverRegistrator< devel::splice::SpliceCreator > reg_SpliceCreator;
static protocols::moves::MoverRegistrator< devel::splice::RBOutMoverCreator > reg_RBOutMoverCreator;
static protocols::moves::MoverRegistrator< devel::splice::RBInMoverCreator > reg_RBInMoverCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::splice::DesignInterfacesOperationCreator > reg_DesignInterfacesOperationCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::splice::FindEndpointsOperationCreator > reg_FindEndpointsOperationCreator;
static protocols::moves::MoverRegistrator< devel::splice::AlignEndsMoverCreator > reg_AlignEndsMoverCreator;
static protocols::moves::MoverRegistrator< devel::cutoutdomain::CutOutDomainCreator > reg_CutOutDomainCreator;

// Task creators
static core::pack::task::operation::TaskOperationRegistrator< devel::denovo_design::task_operations::DesignByResidueCentralityOperationCreator > reg_DesignByResidueCentralityOperationCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::denovo_design::task_operations::DesignCatalyticResiduesOperationCreator > reg_DesignCatalyticResiduesOperationCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::denovo_design::task_operations::DesignByCavityProximityOperationCreator > reg_DesignByCavityProximityOperationCreator;
static core::pack::task::operation::TaskOperationRegistrator< devel::denovo_design::task_operations::DesignRandomRegionOperationCreator > reg_DesignRandomRegionOperationCreator;
core::pack::task::operation::TaskOperationRegistrator< devel::znhash::DisableZnCoordinationResiduesTaskOpCreator > reg_DisableZnCoordinationResiduesTaskOpCreator;

// Filter creators
static protocols::filters::FilterRegistrator< denovo_design::filters::CavityVolumeFilterCreator > reg_CavityVolumeFilterCreator;
static protocols::filters::FilterRegistrator< denovo_design::filters::CoreResiduesPerElementFilterCreator > reg_CoreResiduesPerElementFilterCreator;
static protocols::filters::FilterRegistrator< denovo_design::filters::FoldabilityFilterCreator > reg_FoldabilityFilterCreator;
static protocols::filters::FilterRegistrator< denovo_design::filters::SSShapeComplementarityFilterCreator > reg_SSShapeComplementarityFilterCreator;
static protocols::filters::FilterRegistrator< devel::replica_docking::InteractionScoreFilterCreator > IscCreator_registrator;
static protocols::filters::FilterRegistrator< devel::replica_docking::IrmsdFilterCreator > IrmsdCreator_registrator;
static protocols::filters::FilterRegistrator< devel::replica_docking::FnatFilterCreator > FnatCreator_registrator;
static protocols::filters::FilterRegistrator< devel::replica_docking::LrmsdFilterCreator > LrmsdCreator_registrator;
static protocols::filters::FilterRegistrator< devel::replica_docking::FnonnatFilterCreator > FnonnatCreator_registrator;
static protocols::filters::FilterRegistrator< devel::replica_docking::CaIrmsdFilterCreator > CaIrmsdCreator_registrator;
static protocols::filters::FilterRegistrator< devel::buns::BuriedUnsatHbondFilter2Creator > BuriedUnsatHbondFilter2_registrator;

static protocols::evaluation::EvaluatorRegistrator< devel::replica_docking::WrapFilterAsEvaluatorCreator > reg_WrapFilterAsEvaluatorCreator;

// Energy methods
static core::scoring::methods::EnergyMethodRegistrator< core::scoring::methods::SideChainNeighborsEnergyCreator > reg_SideChainNeighborsEnergyCreator;

void init( int argc, char * argv [] )
{
	register_version_with_core();
	protocols::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	register_version_with_core();
	protocols::init::init( args );
} // init

} // devel

