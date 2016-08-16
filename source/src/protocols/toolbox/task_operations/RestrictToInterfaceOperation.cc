// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceOperation.cc
/// @brief  TaskOperation class that finds an interface and leaves it mobile in the PackerTask
/// @author Steven Lewis smlewi@gmail.com

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <set>

#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToInterfaceOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

/// @details this ctor assumes a pregenerated calculator - if you want a particular non-default cutoff distance
RestrictToInterfaceOperation::RestrictToInterfaceOperation( std::string const & calculator )
: parent(), calculator_name_(calculator)
{
	//I suppose you could reasonably create this object BEFORE the calculator was generated/registered
	//  if( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calculator_name_ ) ){
	//   utility_exit_with_message("In RestrictToInterfaceOperation, calculator " + calculator + " does not exist.");
	//  }
}

/// @brief this ctor will generate the calculator for you (may use defaults)
RestrictToInterfaceOperation::RestrictToInterfaceOperation( core::Size upper_chain, core::Size lower_chain )
: parent(), calculator_name_("")
{
	make_calculator( upper_chain, lower_chain );
}

/// @details private helper function to make calculator - runs in the ctor
void RestrictToInterfaceOperation::make_calculator( core::Size upper_chain, core::Size lower_chain ) {
	make_name( upper_chain, lower_chain );

	using namespace core::pose::metrics;
	if ( CalculatorFactory::Instance().check_calculator_exists( calculator_name_ ) ) {
		Warning() << "In RestrictToInterfaceOperation, calculator " << calculator_name_
			<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		using core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator;
		CalculatorFactory::Instance().register_calculator( calculator_name_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( upper_chain, lower_chain ) ) );
	}
}

/// @details private helper function to name calculator- runs in the ctor
void RestrictToInterfaceOperation::make_name( core::Size upper_chain, core::Size lower_chain ) {
	calculator_name_ = "RTIO_interface_calculator_" + utility::to_string( upper_chain )
		+ '_' + utility::to_string( lower_chain );
}

RestrictToInterfaceOperation::~RestrictToInterfaceOperation() {}

/// @details be warned if you use clone that you'll not get a new interface calculator
core::pack::task::operation::TaskOperationOP RestrictToInterfaceOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterfaceOperation( *this ) );
}

void
RestrictToInterfaceOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	//vector for filling packertask
	utility::vector1_bool repack(pose.total_residue(), false);

	run_calculator(pose, calculator_name_, "interface_residues", repack);

	task.restrict_to_residues(repack);
}

void RestrictToInterfaceOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

core::pack::task::operation::TaskOperationOP
RestrictToInterfaceOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterfaceOperation );
}

void RestrictToInterfaceOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RestrictToInterfaceOperation::provide_xml_schema( xsd );
}

std::string RestrictToInterfaceOperationCreator::keyname() const
{
	return RestrictToInterfaceOperation::keyname();
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
