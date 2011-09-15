// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.cc
/// @brief  TaskOperation class that finds an interface based on:
/// core/pack/task/operation/util/interface_vector_calculate.hh
/// and leaves it mobile in the PackerTask
/// @author Ben Stranges (stranges@unc.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperationCreator.hh>
//#include <protocols/toolbox/task_operations/InterfaceTaskOperation.hh> //shouldn't need this
// Project Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/util/interface_vector_calculate.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/PackerTask.hh>
//#include <protocols/toolbox/pose_metric_calculators/InterfaceVectorDefinitionCalculator.hh>
//#include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <utility/vector1.hh>

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
	parent(), //inits a jump number 1
	jump_active_(true),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{ }

// ///@details this ctor assumes a pregenerated calculator - if you want a particular non-default cutoff distance
// RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( std::string const & calculator )
// 	: parent(), calculator_name_(calculator)
// {}

///@brief this ctor will generate the calculator for you (may use defaults)
///if you want to use chain characters make the calculator that way and pass it to the constructor above
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( core::Size const lower_chain, core::Size const upper_chain ):
	parent(),
	lower_chain_(lower_chain),
	upper_chain_(upper_chain),
	jump_active_(false),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{ }

//full constructor
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( core::Size const lower_chain,
																																				core::Size const upper_chain,
																																				core::Real CB_dist_cutoff,
																																				core::Real nearby_atom_cutoff,
																																				core::Real vector_angle_cutoff,
																																				core::Real vector_dist_cutoff):
	parent(),
	lower_chain_(lower_chain),
	upper_chain_(upper_chain),
	jump_active_(false),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{ }

///@brief this ctor will generate the calculator for you (may use defaults)
///if you want to use chain characters make the calculator that way and pass it to the constructor above
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( utility::vector1_int const movable_jumps ):
	parent(),
	jump_active_(true),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{ set_movable_jumps( movable_jumps ); }

//full constructor
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( utility::vector1_int const movable_jumps ,
																																				core::Real CB_dist_cutoff,
																																				core::Real nearby_atom_cutoff,
																																				core::Real vector_angle_cutoff,
																																				core::Real vector_dist_cutoff):
	parent(),
	jump_active_(true),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{	set_movable_jumps( movable_jumps );}



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


///@details apply function, uses inherited functionality
void
RestrictToInterfaceVectorOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	//using namespace core::pack::task::operation::util;

	//if the jump constructor then itterate through jumps, make union
	if(jump_active_){
		utility::vector1_bool repack_full(pose.total_residue(), false);
		for( utility::vector1_int::const_iterator jj = movable_jumps().begin() ; jj != movable_jumps().end() ; ++jj ){ //itterator fails here
			//std::cout << "Calculating interface for jump: " <<*jj << std::endl;
			//run detection based on jump
			utility::vector1_bool repack =
				core::pack::task::operation::util::calc_interface_vector( pose, *jj,
																																	CB_dist_cutoff_,
																																	nearby_atom_cutoff_,
																																	vector_angle_cutoff_,
																																	vector_dist_cutoff_ );
			//add repack true setting to repack_full
			for(core::Size ii = 1; ii <= repack.size(); ++ii){
				if(repack[ii])
					repack_full[ii] = true;
			}//end add repack to full_repack
		}//itterate over jumps
		task.restrict_to_residues(repack_full);
	}//end if jump active

 	else{ // if using only the two chain case
		//vector for filling packertask
		utility::vector1_bool repack =
			core::pack::task::operation::util::calc_interface_vector( pose,
																																 lower_chain_, upper_chain_,
																																 CB_dist_cutoff_,
																																 nearby_atom_cutoff_,
																																 vector_angle_cutoff_,
																																 vector_dist_cutoff_ );
		task.restrict_to_residues(repack);
	}


}

///@details setters: only exist to pass info from the parser
void
RestrictToInterfaceVectorOperation::upper_chain( core::Size upper_chain){
	upper_chain_ = upper_chain;
	jump_active_ = false;
	//make_name();
}
void
RestrictToInterfaceVectorOperation::lower_chain( core::Size lower_chain){
	lower_chain_ = lower_chain;
	jump_active_ = false;
	//make_name();
}

// void
// RestrictToInterfaceVectorOperation::jump_num( int jump_num ){
// 	add_movable_jump( jump_num );
// 	jump_active_ = true;
// }
//if you want this function use parent class function:
//add_movable_jump( int const additional_jump );

void
RestrictToInterfaceVectorOperation::CB_dist_cutoff( core::Real CB_dist_cutoff){
	CB_dist_cutoff_ = CB_dist_cutoff;
}
void
RestrictToInterfaceVectorOperation::nearby_atom_cutoff(core::Real nearby_atom_cutoff){
	nearby_atom_cutoff_ = nearby_atom_cutoff;
}
void
RestrictToInterfaceVectorOperation::vector_angle_cutoff(core::Real vector_angle_cutoff){
	vector_angle_cutoff_ = vector_angle_cutoff;
}
void
RestrictToInterfaceVectorOperation::vector_dist_cutoff(core::Real vector_dist_cutoff){
	vector_dist_cutoff_ = vector_dist_cutoff;
}

/*
void
RestrictToInterfaceVectorOperation::setup_interface_chains_from_jumps( core::pose::Pose const & pose ){

    // This class wouldn't know what to do with more than jump, so we'll assert that there's only one.
    assert(movable_jumps().size() == 1);
    lower_chain( pose.chain( pose.fold_tree().jump_edge( *it ).start() ) );
    upper_chain( pose.chain( pose.fold_tree().jump_edge( *it ).stop() ) );
}
*/

///@details parse_tag function for the parser, sets reasonable values for all options
void
RestrictToInterfaceVectorOperation::parse_tag( TagPtr tag )
{
	if( ( tag->hasOption("chain1_num" ) || tag->hasOption("chain2_num" ) ) && tag->hasOption("jump" ) )
		utility_exit_with_message( "You can't define chains and jumps" );
	//use chain numbers if given, otherwise use jump
	if( tag->hasOption("chain1_num" ) && tag->hasOption("chain2_num" ) ){
		lower_chain( tag->getOption< core::Size >( "chain1_num") );
		upper_chain( tag->getOption< core::Size >( "chain2_num") );
		jump_active_ = false;
	}
	//go through the jumps, should probably be a utility function for this
	else if( tag->hasOption("jump") )
		{
			//jump_num( tag->getOption< int >( "jump", 1) );
			//get a string of comma separated jumps
			utility::vector1_int jump_vector(0); //init to zero to cause insant problem if nothing passed
			std::vector<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
			for( std::vector<std::string>::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
				// convert to C string, then convert to integer then push back in vector
				int this_jump(std::atoi( it->c_str() ) );
				jump_vector.push_back( this_jump );
				//add_movable_jump (this_jump ); //from parent class
				TR << 	"Adding jump: " << this_jump << ", ";
			}
			TR << std::endl;
			assert( jump_vector.size() > 0 );
			//now set jumps for parent class and locally
			set_movable_jumps( jump_vector );
			jump_active_ = true;
		}//end jumps
	else
		utility_exit_with_message( "Need to define either jump OR (chain1_num AND chain2_num)" );
	//get other possible tags, set default values to something reasonable.
	CB_dist_cutoff( tag->getOption< core::Real >( "CB_dist_cutoff", 10.0) );
	nearby_atom_cutoff( tag->getOption< core::Real >( "nearby_atom_cutoff", 5.5) );
	vector_angle_cutoff( tag->getOption< core::Real >( "vector_angle_cutoff", 75.0) );
	vector_dist_cutoff( tag->getOption< core::Real >( "vector_dist_cutoff", 9.0 ) );

}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
