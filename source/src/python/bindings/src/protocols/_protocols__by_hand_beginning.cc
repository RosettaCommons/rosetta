// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//
///
/// @author Sergey Lyskov
///

#include "boost/python.hpp"

#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/stepwise/options/StepWiseBasicOptions.hh>

#include <utility/factory/WidgetFactory.hh>


// Includes for dummy bindings to simplify import orders
// #include <protocols/jd2/JobOutputterFactory.hh>
// #include <protocols/jd2/JobInputterFactory.hh>
// #include <protocols/jumping/DisulfPairingLibrary.hh>
// #include <protocols/jumping/PairingLibrary.hh>
// #include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.hh>
// #include <protocols/loops/LoopMoverFactory.hh>
// #include <protocols/scoring/InterchainPotential.hh>

void __protocols_by_hand_beginning__()
{
//boost::python::class_< utility::factory::WidgetFactory< protocols::genetic_algorithm::EntityElementCreator >, utility::pointer::owning_ptr< utility::factory::WidgetFactory< protocols::genetic_algorithm::EntityElementCreator > >, boost::noncopyable >( "Py_WidgetFactory_EntityElementCreator", boost::python::no_init);
	boost::python::class_< utility::factory::WidgetFactory< protocols::genetic_algorithm::EntityElementCreator >, boost::noncopyable >( "Py_WidgetFactory_EntityElementCreator", boost::python::no_init );

	//boost::python::class_< protocols::loop_modeling::loggers::Logger, boost::noncopyable >( "__protocols_loop_modeling_loggers_Logger", boost::python::no_init );

    // Dummy bindings to simplify import orders
	boost::python::class_< protocols::loop_modeling::LoopMover, boost::noncopyable >( "__protocols_loop_modeling_LoopMover");

	boost::python::class_< protocols::kinematic_closure::perturbers::Perturber, boost::noncopyable >( "__protocols_kinematic_closure_perturbers_Perturber", boost::python::no_init);

    typedef boost::python::class_< ::protocols::stepwise::options::StepWiseBasicOptions, boost::python::bases< ::basic::resource_manager::ResourceOptions >, boost::noncopyable > StepWiseBasicOptions_exposer_type;
    StepWiseBasicOptions_exposer_type StepWiseBasicOptions_exposer("StepWiseBasicOptions", "protocols/stepwise/options/StepWiseBasicOptions.hh:33", boost::python::init <  >() );

	// Dummy bindings to simplify import orders
	// boost::python::class_< utility::SingletonBase<protocols::jd2::JobOutputterFactory>, boost::noncopyable >( "__utility_SingletonBase_protocols_jd2_JobOutputterFactory__");
	// boost::python::class_< utility::SingletonBase<protocols::jd2::JobInputterFactory>, boost::noncopyable >( "__utility_SingletonBase_protocols_jd2_JobInputterFactory__");
	// boost::python::class_< utility::SingletonBase<protocols::jumping::StandardDisulfPairingLibrary>, boost::noncopyable >( "__utility_SingletonBase_protocols_jumping_StandardDisulfPairingLibrary__");
	// boost::python::class_< utility::SingletonBase<protocols::jumping::StandardPairingLibrary>, boost::noncopyable >( "__utility_SingletonBase_protocols_jumping_StandardPairingLibrary__");
	// boost::python::class_< utility::SingletonBase<protocols::loops::loop_mover::refine::LoopRefineInnerCycleFactory>, boost::noncopyable >( "__utility_SingletonBase_protocols_loops_loop_mover_refine_LoopRefineInnerCycleFactory__");
	// boost::python::class_< utility::SingletonBase<protocols::loops::LoopMoverFactory>, boost::noncopyable >( "__utility_SingletonBase_protocols_loops_LoopMoverFactory__");
	// boost::python::class_< utility::SingletonBase<protocols::scoring::InterchainPotential>, boost::noncopyable >( "__utility_SingletonBase_protocols_scoring_InterchainPotential__");
}
