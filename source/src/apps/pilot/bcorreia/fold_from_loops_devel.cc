// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief


// libRosetta headers

#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/Protocol.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//constraints
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


#include <apps/pilot/bcorreia/fold_from_loops.hh>
#include <devel/fold_from_loops/FoldFromLoopsMover.hh>
#include <devel/fold_from_loops/FoldFromLoops_functions.hh>

// C++ headers
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

// Rosetta Options
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/util.hh>


using namespace core;
using namespace protocols;
using namespace devel::fold_from_loops;
using namespace kinematics;
using namespace basic::options;
using namespace basic::options::OptionKeys;


basic::Tracer TR("bcorreia_fold_from_loops");




// Verifies if it is within the loop or a neighbor range

bool is_loop_neighbor( protocols::loops::Loops & loops, Size & residue, Size & range);


// generates fold trees for topologies with several loops

void fold_tree_generator( protocols::loops::Loops & loops , std::vector<Size> & cutpoints, core::pose::Pose & pose, kinematics::FoldTree & f);



void fold_tree_generator(
		protocols::loops::Loops & loops ,
		std::vector<Size> & cutpoints,
		core::pose::Pose & pose,
		kinematics::FoldTree & f
){
	f.add_edge( 1, pose.total_residue(), Edge::PEPTIDE );

	for (Size i=1;  i < loops.size() ; ++i){
		TR<<"LOOP "<<  i <<std::endl;
		TR<<"Cut point"<< cutpoints[i-1]<<std::endl;
		Size cutpoint= cutpoints[i-1];
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_LOWER, cutpoint); // residue on the pose has to be assigned as a cut
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_UPPER, cutpoint+1);
		Size loop1_midpoint = ((loops[1].stop()-loops[1].start())/2) + loops[1].start();
		TR<<"loop 1 mid point"<< loop1_midpoint<<std::endl;
		Size variable_midpoint = ((loops[i+1].stop()-loops[i+1].start())/2) + loops[i+1].start();
		TR<<"Variable mid_point"<< variable_midpoint <<std::endl;
		f.new_jump( loop1_midpoint, variable_midpoint, cutpoint );
	}

}






////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char* argv [] )
{
	try {

	protocols::abinitio::ClassicAbinitio::register_options();
	protocols::abinitio::AbrelaxApplication::register_options();
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::fragment;
	using namespace protocols;
	using namespace constraints;
	using protocols::jd2::JobDistributor;



	pose::Pose nat_pose;

	pose::Pose extended_pose;


	core::import_pose::pose_from_pdb( nat_pose, basic::options::start_file() );

	extended_pose = nat_pose; //making working copy



	protocols::checkpoint::CheckPointer sliding_checkpoint("closing");


	//Reading loops
	protocols::loops::Loops lr_loops_in;
	if ( option[ OptionKeys::loops::loop_file ].user() ){
		std::string filename( protocols::loops::get_loop_file_name() );

		lr_loops_in.read_loop_file( filename );
	}



	//Defining a fold tree
	kinematics::FoldTree f;
	f.clear();

	std::vector<Size> cut_points;




	//put an option here to protect the input

    core::pose::Pose target_loops;

	std::string swap_loops = option [OptionKeys::fold_from_loops::swap_loops ]().name();

	core::import_pose::pose_from_pdb( target_loops , swap_loops );




	fold_tree_cutpoints_generator(lr_loops_in, cut_points, extended_pose, f);


	ConstraintSetOP ca_cst ( new ConstraintSet() );

	if (option [OptionKeys::fold_from_loops::native_ca_cst ].user() ){

			CA_cst_generator( nat_pose, ca_cst , lr_loops_in, cut_points );


			}





	FragSetOP fragset_large_;
	FragSetOP fragset_small_;

	get_fragments(fragset_large_, fragset_small_ ); // reading_fragments



	FoldFromLoopsMoverOP trial_mover = new FoldFromLoopsMover;

	trial_mover->loops( lr_loops_in );
	trial_mover->set_frag_3mers( fragset_small_ );
	trial_mover->set_frag_9mers( fragset_large_ );
	trial_mover->loops_pdb( target_loops );
	trial_mover->set_ca_csts( ca_cst );



	//add the set up of the



		JobDistributor::get_instance()->go( trial_mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

