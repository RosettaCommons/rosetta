// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--app_name--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#ifdef USEMPI
#include <mpi.h>
#endif

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/moves/mover_schemas.hh>
#include <protocols/jd2/parser/ScoreFunctionLoader.hh>
#include <protocols/jd2/parser/TaskOperationLoader.hh>

#include <protocols/jd3/Job.hh>
#include <protocols/jd3/job_distributors/MPIWorkPartitionJobDistributor.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>

--mover_path--

// utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/datacache/ConstDataMap.hh>
#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


static basic::Tracer TR("--app_name--");

--new_app_options_out--

/// Define your JobQueen
class --class--JobQueen : public protocols::jd3::standard::StandardJobQueen
{
public:
	--class--JobQueen()
	{

		//This is where you will list the options you will be using.  These options are command-line options
		// that can be specified at the <common> and <job> level.

		utility::options::OptionKeyList opts;

		/// Static Functions setup the OptionKeyList:
		//
		//core::scoring::list_read_options_in_get_score_function( opts );
		//core::pack::task::PackerTask::list_options_read( opts );
		//core::pack::task::operation::ReadResfile::list_options_read( opts );
		//protocols::simple_moves::PackRotamersMover::list_options_read( opts );

		add_options( opts );

		/// Options defined in this app:
		//
		//add_option( basic::options::OptionKeys::minimize_sidechains );
		//add_option( basic::options::OptionKeys::min_pack );
		//add_option( basic::options::OptionKeys::off_rotamer_pack );
	}

	~--class--JobQueen() {}

	virtual
	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const
	{
		using namespace utility::tag;
		using namespace protocols::jd2::parser;

		///Here, you specify which tags that are defined in the <common> element.

		/*
		ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
		TaskOperationLoader::provide_xml_schema( job_definition_xsd );

		XMLSchemaSimpleSubelementList subelements;
		subelements
			.add_already_defined_subelement( ScoreFunctionLoader::loader_name(),     & ScoreFunctionLoader::score_function_loader_ct_namer )
			.add_already_defined_subelement( TaskOperationLoader::loader_name(),     & TaskOperationLoader::task_op_loader_ct_namer );
		ct_gen.add_ordered_subelement_set_as_repeatable( subelements );
		*/

	}

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

		/// Here you can specify Schema's to use when creating your XSD for a specific JOB.
		//
		// Example for FixBB JD3 is below:

		/*
		PackRotamersMover::provide_xml_schema(   job_definition_xsd );
		MinMover::provide_xml_schema(            job_definition_xsd );
		ScoreFunctionLoader::provide_xml_schema( job_definition_xsd );
		TaskOperationLoader::provide_xml_schema( job_definition_xsd );

		//task operations and sfxns -- as many as you want
		XMLSchemaSimpleSubelementList task_and_sfxn_subelements;
		task_and_sfxn_subelements
			.add_already_defined_subelement( ScoreFunctionLoader::loader_name(), & ScoreFunctionLoader::score_function_loader_ct_namer )
			.add_already_defined_subelement( TaskOperationLoader::loader_name(), & TaskOperationLoader::task_op_loader_ct_namer );
		job_ct_gen.add_ordered_subelement_set_as_repeatable( task_and_sfxn_subelements );

		//pack -- at most one
		XMLSchemaSimpleSubelementList pack_subelement;
		pack_subelement.add_already_defined_subelement( PackRotamersMoverCreator::mover_name(), & complex_type_name_for_mover );
		job_ct_gen.add_ordered_subelement_set_as_optional( pack_subelement );

		//min -- at most one
		XMLSchemaSimpleSubelementList min_subelement;
		min_subelement.add_already_defined_subelement( MinMoverCreator::mover_name(), & complex_type_name_for_mover );
		job_ct_gen.add_ordered_subelement_set_as_optional( min_subelement );

		*/

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

		/// This is the core of your code.  Most likely, you will need the lines below.
		MoverAndPoseJobOP mature_job( new MoverAndPoseJob );
		core::pose::PoseOP pose = pose_for_job( larval_job, *job_options );
		mature_job->pose( pose );

		TagCOP job_tag;
		if ( larval_job->inner_job()->jobdef_tag() )  {
			job_tag = larval_job->inner_job()->jobdef_tag();
		}

		//Here, you create your mover(s) and call what you need to call. You need to give your protocol a LOCAL
		// Options collection from which to read from.

		--mover_namespace--::--class--OP mover_protocol( new --mover_namespace--::--class--() );
		mover_protocol->initialize_from_options(job_options); //You will need to write this in your protocol!

		// Setup your mover in whatever way you need.
		// SEE the FixBB JD3 example in src/apps/pilot/andrew/fixbb_jd3.cc
		// for more complex setup options, such as combining XML and command-line syntax for specifying job-level options.

		mature_job->mover( mover_protocol );
		return mature_job;
	}

	virtual bool more_jobs_remain() { return false; }


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		--new_app_options_in--
		protocols::jd3::JobDistributorOP jd = protocols::jd3::JobDistributorFactory::create_job_distributor();

		protocols::jd3::JobQueenOP queen( new --class--JobQueen );
		jd->go( queen );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
