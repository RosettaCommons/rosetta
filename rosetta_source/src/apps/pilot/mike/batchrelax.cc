// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/frag_picker/VallChunk.hh>
// AUTO-REMOVED #include <protocols/frag_picker/VallProvider.hh>

#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <boost/cstdint.hpp>
// AUTO-REMOVED #include <boost/unordered_map.hpp>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/init.hh>
// AUTO-REMOVED #include <numeric/HomogeneousTransform.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.hh>
#include <protocols/relax/FastRelax.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>
#include <protocols/match/Hit.fwd.hh>
// AUTO-REMOVED #include <protocols/match/Hit.hh>
// AUTO-REMOVED #include <protocols/match/SixDHasher.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/topology_broker/TopologyBroker.hh>
// AUTO-REMOVED #include <protocols/topology_broker/util.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/fixedsizearray1.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>

// option key includes
// AUTO-REMOVED #include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/mike.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>



static basic::Tracer TR("main");

using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::match;
using namespace protocols::frag_picker;
using core::io::silent::SilentStructFactory;
using core::io::silent::SilentStructOP;




int
main( int argc, char * argv [] )
{
 	using namespace protocols;
 	using namespace protocols::jd2;
 	using namespace basic::options;
 	using namespace basic::options::OptionKeys;
 	using namespace core;
 	using io::silent::SilentStructFactory;
 	using io::silent::SilentStructOP;


 	// initialize core
 	protocols::init(argc, argv);

 	core::Size nstruct = option[ OptionKeys::out::nstruct ];
 	core::Size batch_size = option[ OptionKeys::batch_relax::batch_size ];
 	TR << "The BATCHSIZE: " << batch_size << std::endl;
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

 	io::silent::SilentFileData sfd;
 	std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
 	core::pose::PoseOP native_pose;
 	if( option[ in::file::native ].user() ){
 		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
 	}

 	while( input.has_another_pose() )
 	{
 		// make a list of simple pointers. This should be safe since the input_structs will all remain in scope.
 		core::chemical::ResidueTypeSetCAP rsd_set;
 		if ( option[ in::file::fullatom ]() ) {
 			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
 		} else {
 			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
 		}
 		core::Size count = 0;

 		std::vector < SilentStructOP > input_structs;
 		while( input.has_another_pose() && (count < batch_size ) ) {
 			core::pose::Pose pose;
 			input.fill_pose( pose, *rsd_set );

 			SilentStructOP new_struct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
 			new_struct->fill_struct( pose );
 			input_structs.push_back( new_struct );
			count++;
		}

 		core::scoring::ScoreFunctionOP scorefxn;
 		scorefxn = core::scoring::getScoreFunction();

 		for ( core::Size j=0;j<nstruct;j++ )
 		{
 			std::vector < SilentStructOP > relax_structs;
 			// Make a deep copy of the input_structs list
 			for( std::vector < SilentStructOP >::const_iterator it = input_structs.begin();
 					it != input_structs.end();
 					++ it )
 			{
 				SilentStructOP new_struct;
 				new_struct = (*it)->clone();
 				relax_structs.push_back( new_struct );
 			}

 			protocols::relax::FastRelax relax( scorefxn,  option[ OptionKeys::relax::sequence_file ]() );
 			TR << "BATCHSIZE: " <<  relax_structs.size() << std::endl;
			relax.batch_apply( relax_structs );

 			// Now save the resulting decoys

 			for( std::vector < SilentStructOP >::const_iterator it = relax_structs.begin();
 					it != relax_structs.end();
 					++ it )
 			{
 				if( native_pose ){
 					core::pose::Pose cpose;
 					input.fill_pose( cpose, *rsd_set );
 					core::Real rms = scoring::CA_rmsd( *native_pose, cpose );
 					(*it)->add_energy( "rms", rms, 1.0 );
 				}
 				sfd.write_silent_struct( *(*it), silent_file_ );
 			}
 		} // nstruct for
 	} // while

 	return 0;
}

