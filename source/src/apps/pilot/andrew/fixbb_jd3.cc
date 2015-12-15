// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/public/design/fixbb_jd3.cc
/// @brief  Fixed backbone design, using the jd3 JobDistributor

//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>

#include <core/kinematics/MoveMap.hh>

//protocols library (Movers)
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/MoverAndPoseJob.hh>
#include <protocols/jd3/StandardJobQueen.hh>
#include <protocols/jd3/PoseOutputter.hh>

#include <devel/init.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


class FixbbJobQueen : public protocols::jd3::StandardJobQueen
{
public:
	FixbbJobQueen() {}
	~FixbbJobQueen() {}

	virtual
	protocols::jd3::JobOP
	mature_larval_job( protocols::jd3::LarvalJobCOP job ) {
		using namespace protocols::jd3;
		MoverAndPoseJobOP mature_job( new MoverAndPoseJob );
		core::pose::PoseOP pose = pose_for_job( job );
		mature_job->pose( pose );

		//create a task factory: this will create a new PackerTask for each input pose
		using core::pack::task::operation::TaskOperationCOP;
		core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );

		//create a ScoreFunction from commandline options
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		//create the PackRotamersMover which will do the packing
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );

		pack_mover->task_factory( main_task_factory );
		pack_mover->score_function( score_fxn );

		mature_job->mover( pack_mover );
		return mature_job;
	}

	virtual bool has_job_completed( protocols::jd3::LarvalJobCOP job ) { return pose_outputter().job_has_already_completed( *job ); }
	virtual void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP /*job*/ ) {/*TEMP*/}

	virtual void note_job_completed( protocols::jd3::LarvalJobCOP /*job*/, protocols::jd3::JobStatus /*status*/ ) {}

	virtual void completed_job_result( protocols::jd3::LarvalJobCOP job, protocols::jd3::JobResultOP result ) {
		using namespace protocols::jd3;
		PoseJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( result );
		core::pose::PoseOP pose = pose_result->pose();
		pose_outputter().write_output_pose( *job, *pose );
	}

	virtual
	protocols::jd3::LarvalJobs determine_job_list() {
		protocols::jd3::InnerLarvalJobs inner_job_list = prepare_preliminary_job_list();
		return expand_job_list( inner_job_list );
	}

	virtual bool more_jobs_remain() { return false; }
	virtual std::string job_definition_xsd() const { return ""; }

};


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();
		jd->go( protocols::jd3::JobQueenOP( new FixbbJobQueen ) );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
