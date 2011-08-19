// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.cc
/// @brief  TaskOperation class that finds an interface based on InterfaceVectorDefinitionCalculator and leaves it mobile in the PackerTask
/// @author Ben Stranges (stranges@unc.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/pose_metric_calculators/InterfaceVectorDefinitionCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>


// C++ Headers
#include <set>

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToInterfaceVectorOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

///@details, empty contructor for parser
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation() :
	parent(),
	lower_chain_(1),
	upper_chain_(2),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
	//char_constructor_(false);
{ make_name();}



///@details this ctor assumes a pregenerated calculator - if you want a particular non-default cutoff distance
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( std::string const & calculator )
	: parent(), calculator_name_(calculator)
{}

///@brief this ctor will generate the calculator for you (may use defaults)
///if you want to use chain characters make the calculator that way and pass it to the constructor above
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( core::Size const lower_chain, core::Size const upper_chain ):
	parent(),
	lower_chain_(lower_chain),
	upper_chain_(upper_chain),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{ make_name(); }

//full constructor
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( core::Size const lower_chain,
																																				core::Size upper_chain,
																																				core::Real CB_dist_cutoff,
																																				core::Real nearby_atom_cutoff,
																																				core::Real vector_angle_cutoff,
																																				core::Real vector_dist_cutoff):
	parent(),
	lower_chain_(lower_chain),
	upper_chain_(upper_chain),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{ make_name();}

//class member functions
RestrictToInterfaceVectorOperation::~RestrictToInterfaceVectorOperation() {}

core::pack::task::operation::TaskOperationOP
RestrictToInterfaceVectorOperationCreator::create_task_operation() const
{
	return new RestrictToInterfaceVectorOperation;
}

///@details be warned if you use clone that you'll not get a new interface calculator
core::pack::task::operation::TaskOperationOP RestrictToInterfaceVectorOperation::clone() const
{
	return new RestrictToInterfaceVectorOperation( *this );
}

///@details private helper function to make calculator - runs in the ctor
void RestrictToInterfaceVectorOperation::make_calculator( const std::string & calculator_name ) const {
	using namespace core::pose::metrics;
	using namespace protocols::toolbox::pose_metric_calculators;

	if( CalculatorFactory::Instance().check_calculator_exists( calculator_name ) ){
		Warning() << "In RestrictToInterfaceVectorOperation, calculator " << calculator_name
							<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
		InterfaceVectorDefinitionCalculatorOP interface_calc = new InterfaceVectorDefinitionCalculator( lower_chain_, upper_chain_, CB_dist_cutoff_, nearby_atom_cutoff_, vector_angle_cutoff_, vector_dist_cutoff_ );

		CalculatorFactory::Instance().register_calculator( calculator_name, interface_calc );
	}
}

///@details private helper function to name calculator will be unique unless the setters are altered
void RestrictToInterfaceVectorOperation::make_name() {
	calculator_name_ = "RIVO_interface_calculator_" + utility::to_string( lower_chain_ )
		+ '_' + utility::to_string( upper_chain_ ) +  '_' +
		utility::to_string( CB_dist_cutoff_ ) + '_' +
		utility::to_string( nearby_atom_cutoff_ ) + '_' +
		utility::to_string( vector_angle_cutoff_ ) + '_' +
		utility::to_string( vector_dist_cutoff_ );
}
///@details apply function, uses inherited functionality
void
RestrictToInterfaceVectorOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{

	//can do this here if we want setters
	make_calculator(calculator_name_);

 	//vector for filling packertask
 	utility::vector1_bool repack(pose.total_residue(), false);
	run_calculator(pose, calculator_name_, "interface_residues", repack);
 	task.restrict_to_residues(repack);

}

///@details setters: only exist to pass info to the calculator
/// calculator name must be remade when these are accessed
void
RestrictToInterfaceVectorOperation::upper_chain( core::Size upper_chain){
	upper_chain_ = upper_chain;
	make_name();
}
void
RestrictToInterfaceVectorOperation::lower_chain( core::Size lower_chain){
	lower_chain_ = lower_chain;
	make_name();
}
void
RestrictToInterfaceVectorOperation::CB_dist_cutoff( core::Real CB_dist_cutoff){
	CB_dist_cutoff_ = CB_dist_cutoff;
	make_name();
}
void
RestrictToInterfaceVectorOperation::nearby_atom_cutoff(core::Real nearby_atom_cutoff){
	nearby_atom_cutoff_ = nearby_atom_cutoff;
	make_name();
}
void
RestrictToInterfaceVectorOperation::vector_angle_cutoff(core::Real vector_angle_cutoff){
	vector_angle_cutoff_ = vector_angle_cutoff;
	make_name();
}
void
RestrictToInterfaceVectorOperation::vector_dist_cutoff(core::Real vector_dist_cutoff){
	vector_dist_cutoff_ = vector_dist_cutoff;
	make_name();
}

///@details parse_tag function for the parser, sets reasonable values for all options
void
RestrictToInterfaceVectorOperation::parse_tag( TagPtr tag )
{
	lower_chain( tag->getOption< core::Size >( "chain1_num", 1) );
	upper_chain( tag->getOption< core::Size >( "chain2_num", 2) );
	CB_dist_cutoff( tag->getOption< core::Real >( "CB_dist_cutoff", 10.0) );
	nearby_atom_cutoff( tag->getOption< core::Real >( "nearby_atom_cutoff", 5.5) );
	vector_angle_cutoff( tag->getOption< core::Real >( "vector_angle_cutoff", 75.0) );
	vector_dist_cutoff( tag->getOption< core::Real >( "vector_dist_cutoff", 9.0 ) );
}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
