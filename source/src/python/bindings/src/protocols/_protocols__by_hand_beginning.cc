
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


//#include <protocols/loop_modeling/loggers/Logger.hh>



#include <utility/factory/WidgetFactory.hh>

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
}
