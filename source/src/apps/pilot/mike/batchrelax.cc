// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/frag_picker/VallChunk.hh>

#include <utility/pointer/owning_ptr.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>

#include <devel/init.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/batch_relax.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "main" );

int
main( int argc, char * argv [] )
{
	try {

		using namespace core;
		using namespace protocols;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using io::silent::SilentStructFactory;
		using io::silent::SilentStructOP;


		// initialize core
		devel::init(argc, argv);

		evaluation::PoseEvaluatorsOP evaluators_( new protocols::evaluation::PoseEvaluators() );
		evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluators_);


		core::Size nstruct = option[ OptionKeys::out::nstruct ];
		core::Size batch_size = option[ OptionKeys::batch_relax::batch_size ];
		TR << "The BATCHSIZE: " << batch_size << std::endl;
		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

		io::silent::SilentFileOptions opts; // initialized from the command line
		io::silent::SilentFileData sfd( opts );
		std::string silent_file_ = option[ OptionKeys::out::file::silent ]();

		//core::Size struct_count = 0;  // unused ~Labonte
		while ( input.has_another_pose() )
				{
			// make a list of simple pointers. This should be safe since the input_structs will all remain in scope.
			core::chemical::ResidueTypeSetCOP rsd_set;
			if ( option[ in::file::fullatom ]() ) {
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			} else {
				rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
			}
			core::Size count = 0;

			std::vector < SilentStructOP > input_structs;
			while ( input.has_another_pose() && (count < batch_size ) ) {
				core::pose::Pose pose;
				input.fill_pose( pose, *rsd_set );

				SilentStructOP new_struct = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
				new_struct->fill_struct( pose );
				input_structs.push_back( new_struct );
				count++;
			}

			core::scoring::ScoreFunctionOP scorefxn;
			scorefxn = core::scoring::get_score_function();

			for ( core::Size j=0; j<nstruct; j++ ) {
				std::vector < SilentStructOP > relax_structs;
				// Make a deep copy of the input_structs list
				for ( std::vector < SilentStructOP >::const_iterator it = input_structs.begin();
						it != input_structs.end();
						++ it ) {
					SilentStructOP new_struct;
					new_struct = (*it)->clone();
					relax_structs.push_back( new_struct );
				}

				protocols::relax::FastRelax relax( scorefxn,  option[ OptionKeys::relax::sequence_file ]() );
				TR << "BATCHSIZE: " <<  relax_structs.size() << std::endl;
				long starttime = time(NULL);
				relax.batch_apply( relax_structs );
				long endtime = time(NULL);
				TR << "TIME: " << endtime - starttime << " seconds" << std::endl;

				// Now save the resulting decoys

				for ( std::vector < SilentStructOP >::const_iterator it = relax_structs.begin();
						it != relax_structs.end();
						++ it ) {
					if ( evaluators_->size() ) {
						core::pose::Pose cpose;
						input.fill_pose( cpose, *rsd_set );
						evaluators_->apply( cpose, "tag" , *(*it) );
					}


					sfd.write_silent_struct( *(*it), silent_file_ );
				}
			} // nstruct for
		} // while

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
