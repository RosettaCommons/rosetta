// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /src/apps/public/design/kellogg_ddg_monomer_jd3.cc
/// @brief  The protocol for calculating ddGs of mutation for monomers developed by Liz Kellogg
///         using the jd3 JobDistributor to parallelize the computations

#ifdef USEMPI
#include <mpi.h>
#endif

//core library
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>

#include <core/kinematics/MoveMap.hh>

//protocols library (Movers)
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinPackMoverCreator.hh>
#include <protocols/simple_moves/PackRotamersMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/MinMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/util.hh>

#include <protocols/jd2/parser/ScoreFunctionLoader.hh>
#include <protocols/jd2/parser/TaskOperationLoader.hh>

#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/MoverAndPoseJob.hh>
#include <protocols/jd3/StandardJobQueen.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>

#include <devel/init.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

//local options
namespace basic { namespace options { namespace OptionKeys {
basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
basic::options::BooleanOptionKey const min_pack("min_pack");
basic::options::BooleanOptionKey const off_rotamer_pack("off_rotamer_pack");
}}}//basic::options::OptionKeys


static THREAD_LOCAL basic::Tracer TR( "apps.pilot.andrew.fixbb_jd3" );

class FixbbJobQueen : public protocols::jd3::StandardJobQueen
{
public:
	FixbbJobQueen()
	{
		utility::options::OptionKeyList opts;
		core::scoring::list_read_options_in_get_score_function( opts );
		core::pack::task::PackerTask::list_options_read( opts );
		core::pack::task::operation::ReadResfile::list_options_read( opts );
		protocols::simple_moves::PackRotamersMover::list_options_read( opts );
		add_options( opts );
		add_option( basic::options::OptionKeys::minimize_sidechains );
		add_option( basic::options::OptionKeys::min_pack );
		add_option( basic::options::OptionKeys::off_rotamer_pack );
	}

	~FixbbJobQueen() {}

	virtual
	void append_job_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen
	) const
	{
		using namespace utility::tag;
		using namespace protocols::simple_moves;
		using namespace protocols::moves;
		using namespace protocols::jd2::parser;

		PackRotamersMover::provide_xml_schema(   job_definition_xsd );
		MinMover::provide_xml_schema(            job_definition_xsd );
		ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
		TaskOperationLoader::provide_xml_schema( job_definition_xsd );

		// task operations and sfxns -- as many as you want
		XMLSchemaSimpleSubelementList task_and_sfxn_subelements;
		task_and_sfxn_subelements
			.add_already_defined_subelement( ScoreFunctionLoader::loader_name(), & ScoreFunctionLoader::score_function_loader_ct_namer )
			.add_already_defined_subelement( TaskOperationLoader::loader_name(), & TaskOperationLoader::task_op_loader_ct_namer );
		job_ct_gen.add_ordered_subelement_set_as_repeatable( task_and_sfxn_subelements );

		// pack -- at most one
		XMLSchemaSimpleSubelementList pack_subelement;
		pack_subelement.add_already_defined_subelement( PackRotamersMoverCreator::mover_name(), & complex_type_name_for_mover );
		job_ct_gen.add_ordered_subelement_set_as_optional( pack_subelement );

		// min -- at most one
		XMLSchemaSimpleSubelementList min_subelement;
		min_subelement.add_already_defined_subelement( MinMoverCreator::mover_name(), & complex_type_name_for_mover );
		job_ct_gen.add_ordered_subelement_set_as_optional( min_subelement );
	}

	virtual
	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const
	{
		using namespace utility::tag;
		using namespace protocols::jd2::parser;
		ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
		TaskOperationLoader::provide_xml_schema( job_definition_xsd );

		XMLSchemaSimpleSubelementList subelements;
		subelements
			.add_already_defined_subelement( ScoreFunctionLoader::loader_name(),     & ScoreFunctionLoader::score_function_loader_ct_namer )
			.add_already_defined_subelement( TaskOperationLoader::loader_name(),     & TaskOperationLoader::task_op_loader_ct_namer );
		ct_gen.add_ordered_subelement_set_as_repeatable( subelements );
		// TO DO: Is there anything else that needs to be set?
	}


	virtual
	protocols::jd3::JobOP
	complete_larval_job_maturation(
		protocols::jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< protocols::jd3::JobResultCOP > const &
	)
	{

		TR << "Completing larval job maturation" << std::endl;

		using namespace protocols::jd3;
		using namespace protocols::moves;
		using namespace protocols::simple_moves;
		using namespace utility::tag;
		using namespace basic::datacache;
		using namespace core::pack::task::operation;

		MoverAndPoseJobOP mature_job( new MoverAndPoseJob );
		core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
		mature_job->pose( pose );


		TagCOP job_tag;
		if ( larval_job->inner_job()->jobdef_tag() )  {
			job_tag = larval_job->inner_job()->jobdef_tag();
		}
		SequenceMoverOP seq( new SequenceMover );

		if ( job_tag ) {
			TR << "Initializing fixbb job from job_tag" << std::endl;
			// parse the score functions and task operations in the common & job tags
			DataMap datamap;
			TagCOP common_tag = common_block_tags();
			parse_sfxns_and_taskops( *pose, common_tag, datamap );
			parse_sfxns_and_taskops( *pose,    job_tag, datamap );

			MoverOP pack_mover;
			if ( job_tag->hasTag( PackRotamersMoverCreator::mover_name() ) ) {
				pack_mover = MoverOP( new PackRotamersMover( *job_options ));
				TR << "Calling PRM::parse_my_tag" << std::endl;
				pack_mover->parse_my_tag(
					job_tag->getTag( PackRotamersMoverCreator::mover_name() ),
					datamap,
					Mover::Filters_map(),
					Movers_map(),
					*pose );
			} else {
				using namespace core::pack::task::operation;
				// read the score function from the job options object; create a task operation to initialize
				// the (eventually created) packer task from the job options as well.
				PackRotamersMoverOP prm( new PackRotamersMover( *job_options ));
				pack_mover = prm;
				prm->score_function( core::scoring::get_score_function( *job_options ) );
				core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
				task_factory->push_back( TaskOperationOP( new InitializeFromOptionCollection( job_options )));
				task_factory->push_back( TaskOperationOP( new ReadResfile( *job_options )));
				prm->task_factory( task_factory );
			}
			seq->add_mover( pack_mover );

			MoverOP min_mover;
			if ( job_tag->hasTag( MinMoverCreator::mover_name() ) ) {
				min_mover = MoverOP( new MinMover );
				min_mover->parse_my_tag(
					job_tag->getTag( MinMoverCreator::mover_name()),
					datamap,
					Mover::Filters_map(),
					Movers_map(),
					*pose );
			} else if ( (*job_options)[ basic::options::OptionKeys::minimize_sidechains ] ) {
				core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function( *job_options );
				core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
				movemap->set_chi( true );
				min_mover = MinMoverOP( new MinMover(
					movemap,
					score_fxn,
					(*job_options)[ basic::options::OptionKeys::run::min_type ].value(),
					0.01,
					true
					));
			}
			if ( min_mover ) seq->add_mover( min_mover );

		} else {
			// initialize the pack rotamers mover and (possibly) the min mover from the job-options object.
			// create a task factory and initialize w/ a read-resfile and init from options pair of task operations.
			core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
			main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromOptionCollection( job_options ) ));
			main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );

			//create a ScoreFunction from commandline options
			core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function( *job_options );

			//create the PackRotamersMover which will do the packing
			protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
			pack_mover->task_factory( main_task_factory );
			pack_mover->score_function( score_fxn );
			seq->add_mover( pack_mover );

			MoverOP min_mover;
			if ( (*job_options)[ basic::options::OptionKeys::minimize_sidechains ] ) {
				core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function( *job_options );
				core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
				movemap->set_chi( true );
				min_mover = MinMoverOP( new MinMover(
					movemap,
					score_fxn,
					(*job_options)[ basic::options::OptionKeys::run::min_type ].value(),
					0.01,
					true
					));
			}
			if ( min_mover ) seq->add_mover( min_mover );

		}

		mature_job->mover( seq );
		return mature_job;
	}

	//virtual bool has_job_completed( protocols::jd3::LarvalJobCOP job ) { return pose_outputter_for_job( *job->inner_job() )->job_has_already_completed( *job ); }
	virtual void mark_job_as_having_begun( protocols::jd3::LarvalJobCOP /*job*/ ) {/*TEMP*/}

	//virtual void note_job_completed( protocols::jd3::LarvalJobCOP /*job*/, protocols::jd3::JobStatus /*status*/ ) {}

	//virtual void completed_job_result( protocols::jd3::LarvalJobCOP job, protocols::jd3::JobResultOP result ) {
	// using namespace protocols::jd3;
	// PoseJobResultOP pose_result = utility::pointer::dynamic_pointer_cast< PoseJobResult > ( result );
	// core::pose::PoseOP pose = pose_result->pose();
	// utility::options::OptionCollectionCOP job_options = options_for_job( *job->inner_job() );
	// pose_outputter_for_job( *job->inner_job() )->write_output_pose( *job, *job_options, *pose );
	//}

	virtual bool more_jobs_remain() { return false; }

	void
	parse_sfxns_and_taskops(
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const
	{
		if ( ! tag ) return;

		using namespace utility::tag;
		using namespace protocols::jd2::parser;

		for ( Tag::tags_t::const_iterator iter = tag->getTags().begin(); iter != tag->getTags().end(); ++iter ) {
			if ( (*iter)->getName() == ScoreFunctionLoader::loader_name() ) {
				ScoreFunctionLoader sfxn_loader;
				sfxn_loader.load_data( pose, *iter, datamap );
			} else if ( (*iter)->getName() == TaskOperationLoader::loader_name() ) {
				TaskOperationLoader taskop_loader;
				taskop_loader.load_data( pose, *iter, datamap );
			}
		}

	}

};


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		option.add( minimize_sidechains, "Do minimization of side chains after rotamer packing").def(false);
		option.add( min_pack, "Pack and minimize sidechains simultaneously").def(false);
		option.add( off_rotamer_pack, "Pack using a continuous sidechain rotamer library").def(false);

		devel::init(argc, argv);

		protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();
		protocols::jd3::JobQueenOP queen( new FixbbJobQueen );
		//std::cout << "Fixbb job definition file\n" << queen->job_definition_xsd() << std::endl;
		jd->go( queen );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
