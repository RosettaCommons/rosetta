// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilot/andrew/relax_jd3.cc
/// @brief  Relax, using the jd3 JobDistributor

#ifdef USEMPI
#include <mpi.h>
#endif

// Package Headers
#include <protocols/relax/relax_main.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

//core library
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pack/task/operation/TaskOperation.hh>
//#include <core/pack/task/operation/TaskOperations.hh>
//#include <core/pack/task/TaskFactory.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>

//#include <core/kinematics/MoveMap.hh>

//protocols library (Movers)
//#include <protocols/minimization_packing/MinPackMover.hh>
//#include <protocols/minimization_packing/PackRotamersMover.hh>
//#include <protocols/minimization_packing/MinPackMoverCreator.hh>
//#include <protocols/minimization_packing/PackRotamersMoverCreator.hh>
//#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
//#include <protocols/minimization_packing/MinMover.hh>
//#include <protocols/minimization_packing/MinMoverCreator.hh>
//#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
//#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/mover_schemas.hh>

#include <protocols/jd2/util.hh>
//#include <protocols/jd2/parser/TaskOperationLoader.hh>

#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
//#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>


#include <devel/init.hh>

// Basic headers
//#include <basic/datacache/ConstDataMap.hh>
//#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


// option key includes
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

//local options
//namespace basic { namespace options { namespace OptionKeys {
//basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
//basic::options::BooleanOptionKey const min_pack("min_pack");
//basic::options::BooleanOptionKey const off_rotamer_pack("off_rotamer_pack");
//}}}//basic::options::OptionKeys


static basic::Tracer TR( "apps.pilot.andrew.relax_jd3" );

class RelaxJobQueen : public protocols::jd3::standard::StandardJobQueen
{
public:
	RelaxJobQueen()
	{
		using namespace basic::options::OptionKeys;
		utility::options::OptionKeyList opts;
		core::scoring::list_read_options_in_get_score_function( opts );
		protocols::relax::options_for_generate_relax_from_cmd( opts );
		protocols::jd2::options_for_set_native_in_mover( opts );
		core::scoring::constraints::options_for_get_cst_file_option( opts );
		core::scoring::constraints::options_for_get_cst_fa_file_option( opts );
		protocols::simple_moves::symmetry::SetupForSymmetryMover::options_read_in_ctor( opts );
		protocols::electron_density::SetupForDensityScoringMover::options_read_in_ctor( opts );
		opts
			+ constraints::cst_fa_file
			+ constraints::cst_file
			+ edensity::mapfile
			+ symmetry::symmetry_definition
			+ relax::superimpose_to_file
			+ relax::superimpose_to_native;
		add_options( opts );
	}

	~RelaxJobQueen() override = default;




	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< protocols::jd3::JobResultCOP > const &
	) override
	{

		TR << "Completing larval job maturation" << std::endl;

		using namespace protocols::jd3;
		using namespace protocols::jd3::standard;
		using namespace protocols::moves;
		using namespace protocols::simple_moves;
		using namespace utility::tag;
		using namespace basic::datacache;
		using namespace basic::options;
		using namespace core::pack::task::operation;

		MoverAndPoseJobOP mature_job( new MoverAndPoseJob );
		core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
		mature_job->pose( pose );


		TagCOP job_tag;
		if ( larval_job->inner_job()->jobdef_tag() )  {
			job_tag = larval_job->inner_job()->jobdef_tag();
		}

		//if ( job_tag ) {
		// // there should be something here
		//} else {
		// This code belongs in the else once we allow real data into the
		// job tag, but for now, everything will be loaded through the
		// options system.
		// You might reasonably wonder why this code is a near duplication of the
		// Relax_main function. I guess my thought is that this code will eventually
		// get better by reading data out of tags, and to do that properly, it will
		// need to diverge from Relax_main.

		protocols::moves::MoverOP protocol = protocols::relax::generate_relax_from_cmd( *job_options, false );
		protocols::jd2::set_native_in_mover( *protocol, *job_options );


		// add constraints from cmd line
		if ( (*job_options)[ OptionKeys::constraints::cst_fa_file ].user() ||
				(*job_options)[ OptionKeys::constraints::cst_file ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			protocols::simple_moves::ConstraintSetMoverOP load_csts( new protocols::simple_moves::ConstraintSetMover );
			if ( (*job_options)[ OptionKeys::constraints::cst_fa_file ].user() ) {
				load_csts->constraint_file( core::scoring::constraints::get_cst_fa_file_option( *job_options ) );
			} else {
				load_csts->constraint_file( core::scoring::constraints::get_cst_file_option( *job_options ) );
			}
			seqmov->add_mover( load_csts );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// set pose for density scoring if a map was input
		// (potentially) dock map into density
		if ( (*job_options)[ OptionKeys::edensity::mapfile ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::electron_density::SetupForDensityScoringMover( *job_options ) ) );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// setup symmetry mover ... this should happen _before_ SetupForDensityScoringMover
		//   to avoid adding extra VRTs
		if ( (*job_options)[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::simple_moves::symmetry::SetupForSymmetryMover( *job_options ) ) );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}

		// superimpose input model to the native structure or to a supplied PDB file
		if ( (*job_options)[ OptionKeys::relax::superimpose_to_file ].user() ||
				(*job_options)[ OptionKeys::relax::superimpose_to_native ].user()
				) {
			core::pose::Pose ref_pose;
			std::string ref_filename;
			if (  (*job_options)[ OptionKeys::relax::superimpose_to_file ].user() ) ref_filename = (*job_options)[ basic::options::OptionKeys::relax::superimpose_to_file ]();
			if (  (*job_options)[ OptionKeys::relax::superimpose_to_native ].user() ) ref_filename =  (*job_options)[ basic::options::OptionKeys::in::file::native ]();

			// Arguably there is a missing feature here of letting the derived
			// job queen to have the standard job queen store the native pose after
			// it is first read in. The way this stands, the native pose is loaded
			// over and over.
			core::import_pose::pose_from_file( ref_pose, ref_filename , core::import_pose::PDB_file);
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			protocols::simple_moves::SuperimposeMoverOP sm( new protocols::simple_moves::SuperimposeMover );
			sm->set_reference_pose( ref_pose );
			seqmov->add_mover( sm );
			seqmov->add_mover( protocol );
			protocol = seqmov;
		}



		//} end else-job-tag

		mature_job->mover( protocol );
		return mature_job;
	}

	//virtual bool has_job_completed( protocols::jd3::LarvalJobCOP job ) { return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job ); }
	void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP /*job*/ ) override {/*TEMP*/}

	virtual bool more_jobs_remain() { return false; }

};


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();
		protocols::jd3::JobQueenOP queen( new RelaxJobQueen );
		//std::cout << "Fixbb job definition file\n" << queen->job_definition_xsd() << std::endl;
		jd->go( queen );
	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
