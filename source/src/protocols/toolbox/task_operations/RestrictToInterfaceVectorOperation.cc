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
#include <core/conformation/Conformation.hh>
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
// AUTO-REMOVED #include <set>

#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToInterfaceVectorOperation" );

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
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation( core::Size const lower_chain_id, core::Size const upper_chain_id ):
	parent(),
	jump_active_(false),
	CB_dist_cutoff_( 11.0 ),
	nearby_atom_cutoff_( 5.5 ),
	vector_angle_cutoff_( 75.0 ),
	vector_dist_cutoff_( 9.0 )
{
	lower_chain(lower_chain_id);
	upper_chain(upper_chain_id);

}

//full constructor
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation(
	core::Size const lower_chain_id,
	core::Size const upper_chain_id,
	core::Real CB_dist_cutoff,
	core::Real nearby_atom_cutoff,
	core::Real vector_angle_cutoff,
	core::Real vector_dist_cutoff
):
	parent(),
	jump_active_(false),
	CB_dist_cutoff_( CB_dist_cutoff ),
	nearby_atom_cutoff_( nearby_atom_cutoff ),
	vector_angle_cutoff_( vector_angle_cutoff ),
	vector_dist_cutoff_( vector_dist_cutoff )
{
	lower_chain(lower_chain_id);
	upper_chain(upper_chain_id);
}

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
RestrictToInterfaceVectorOperation::RestrictToInterfaceVectorOperation(
	utility::vector1_int const movable_jumps ,
	core::Real CB_dist_cutoff,
	core::Real nearby_atom_cutoff,
	core::Real vector_angle_cutoff,
	core::Real vector_dist_cutoff
):
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
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterfaceVectorOperation );
}

///@details be warned if you use clone that you'll not get a new interface calculator
core::pack::task::operation::TaskOperationOP RestrictToInterfaceVectorOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new RestrictToInterfaceVectorOperation( *this ) );
}


///@details apply function, uses inherited functionality
void
RestrictToInterfaceVectorOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{

	//if the jump constructor then itterate through jumps, make union
	if(jump_active_){
		utility::vector1_bool repack_full(pose.total_residue(), false);
		for( utility::vector1_int::const_iterator jj = movable_jumps().begin() ; jj != movable_jumps().end() ; ++jj ){
			//std::cout << "Calculating interface for jump: " <<*jj << std::endl;
			//run detection based on jump
			if( *jj == 0 || core::Size(*jj) > pose.num_jump() ) {
				TR.Error << "In RestrictToInterfaceVectorOperation: Jump " << *jj << " not found in pose with " << pose.num_jump() << " jumps." << std::endl;
				utility_exit_with_message("Specified jump does not exist.");
			}
			utility::vector1_bool repack =
				core::pack::task::operation::util::calc_interface_vector(
					pose,
					*jj,
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
 		utility::vector1_bool repack_full(pose.total_residue(),false);

 		for(utility::vector1<core::Size>::const_iterator lower_chain_it = lower_chains_.begin();
				lower_chain_it != lower_chains_.end(); ++lower_chain_it)
			{
				core::Size current_lower_chain = *lower_chain_it;

				for(utility::vector1<core::Size>::const_iterator upper_chain_it = upper_chains_.begin();
						upper_chain_it != upper_chains_.end(); ++upper_chain_it)
					{
						core::Size current_upper_chain = *upper_chain_it;
						TR.Debug << "calculating_interface between: " << current_lower_chain << " " << current_upper_chain <<std::endl;
						if( current_lower_chain == 0 || current_lower_chain > pose.conformation().num_chains() ||
								current_upper_chain == 0 || current_upper_chain > pose.conformation().num_chains() ) {
							TR.Error << "In RestrictToInterfaceVectorOperation: Cannot find interface between chains " << current_lower_chain <<
									" and " << current_upper_chain << ". Pose only has " << pose.conformation().num_chains() << " chains." << std::endl;
							utility_exit_with_message("Specified chain does not exist.");
						}
						utility::vector1_bool repack =
							core::pack::task::operation::util::calc_interface_vector( pose,
								current_lower_chain, current_upper_chain,
								CB_dist_cutoff_,
								nearby_atom_cutoff_,
								vector_angle_cutoff_,
								vector_dist_cutoff_ );
						for(core::Size ii = 1; ii <=repack.size(); ++ii)
							{
								if(repack[ii])
									{
										repack_full[ii] = true;
									}
							}
					}
			}

		task.restrict_to_residues(repack_full);
	}


}

///@details setters: only exist to pass info from the parser
void
RestrictToInterfaceVectorOperation::upper_chain( core::Size upper_chain){
	upper_chains_.clear();
	upper_chains_.push_back(upper_chain);
	jump_active_ = false;
	assert(upper_chains_.size() == 1);
	//make_name();
}

void RestrictToInterfaceVectorOperation::upper_chain(utility::vector1<core::Size> upper_chain)
{
	upper_chains_ = upper_chain;
	jump_active_ = false;
}

void
RestrictToInterfaceVectorOperation::lower_chain( core::Size lower_chain){
	lower_chains_.clear();
	lower_chains_.push_back(lower_chain);
	jump_active_ = false;
	assert(lower_chains_.size() == 1);
	//make_name();
}


void RestrictToInterfaceVectorOperation::lower_chain(utility::vector1<core::Size> lower_chain)
{
	lower_chains_ = lower_chain;
	jump_active_ = false;
}

//if you want to change the jump use parent class function:
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
RestrictToInterfaceVectorOperation::parse_tag( TagCOP tag , DataMap & )
{
	if( ( tag->hasOption("chain1_num" ) || tag->hasOption("chain2_num" ) ) && tag->hasOption("jump" ) )
		utility_exit_with_message( "You can't define chains and jumps" );
	//use chain numbers if given, otherwise use jump
	if( tag->hasOption("chain1_num" ) && tag->hasOption("chain2_num" ) ){

		utility::vector1<core::Size> lower_chains = utility::string_split(tag->getOption<std::string>("chain1_num"),',',core::Size());
		utility::vector1<core::Size> upper_chains = utility::string_split(tag->getOption<std::string>("chain2_num"),',',core::Size());

		lower_chain( lower_chains );
		upper_chain( upper_chains );
		jump_active_ = false;
	}
	//go through the jumps, should probably be a utility function for this
	else if( tag->hasOption("jump") )
		{
			//jump_num( tag->getOption< int >( "jump", 1) );
			//get a string of comma separated jumps
			utility::vector1_int jump_vector(0); //init to zero to cause insant problem if nothing passed
			utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );
			for( utility::vector1<std::string>::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
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
