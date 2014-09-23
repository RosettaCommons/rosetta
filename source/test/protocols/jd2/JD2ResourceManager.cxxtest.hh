// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/jd2/JD2ResourceManager.cxxtest.hh
/// @brief test suite for protocols::jd2::JD2ResourceManager and protocols::resource_manager::LazyResourceManager
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/loops/LoopsFileOptions.hh>

// Package headers
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/FallbackConfigurationCreator.hh>
#include <basic/resource_manager/FallbackConfigurationFactory.hh>
#include <basic/resource_manager/LazyResourceManager.hh>
#include <basic/resource_manager/ResourceManagerFactory.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <basic/resource_manager/ResourceOptionsFactory.hh>
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderCreator.hh>
#include <basic/resource_manager/ResourceLoaderFactory.hh>
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/ResourceLocatorCreator.hh>
#include <basic/resource_manager/ResourceLocatorFactory.hh>
#include <basic/resource_manager/types.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// Numberic headers

// C++ headers
#include <string>

using namespace basic::resource_manager;
using namespace protocols::jd2;

class DummyResource : public utility::pointer::ReferenceCount
{
public:
	int somevar_;
};

class DummyResourceOptions : public ResourceOptions {
public:
	DummyResourceOptions() : somevar_( 1 ) {}
	DummyResourceOptions(int somevar) : somevar_( somevar ) {}

	virtual
	void
	parse_my_tag( utility::tag::TagCOP tag )
	{
		somevar_ = tag->getOption< int >( "somevar", 1 );
	}

	virtual
	std::string
	type() const { return "DummyResourceOptions"; }

public:
	int somevar_;

};

typedef utility::pointer::shared_ptr< DummyResourceOptions > DummyResourceOptionsOP;

class DummyResourceOptionsCreator : public ResourceOptionsCreator {
public:
	virtual std::string options_type() const { return "DummyResourceOptions"; }
	virtual ResourceOptionsOP create_options() const { return ResourceOptionsOP( new DummyResourceOptions ); }
};



class DummyResourceLoader : public ResourceLoader {
public:
	DummyResourceLoader() {}

	virtual
	utility::pointer::ReferenceCountOP
	create_resource(
		ResourceOptions const & opts,
		std::string const & ,
		std::istream &
	) const {
		DummyResourceOptions const * dopts = dynamic_cast< DummyResourceOptions const * > ( &opts );
		if ( ! dopts ) {
			utility_exit_with_message( "Bad unit test" );
		}
		DummyResource * dummy_ptr = new DummyResource;
		dummy_ptr->somevar_ = dopts->somevar_;
		return utility::pointer::ReferenceCountOP(dummy_ptr);
	}

	virtual
	ResourceOptionsOP
	default_options() const {
		return ResourceOptionsOP( new DummyResourceOptions );
	}

};

class DummyResourceLoaderCreator : public ResourceLoaderCreator {
public:
	virtual
	ResourceLoaderOP
	create_resource_loader() const {
		return ResourceLoaderOP( new DummyResourceLoader );
	}

	virtual
	std::string loader_type() const { return "DummyResource"; }

};

class DummyResourceStream : public ResourceStream
{
public:
	virtual
	std::istream &
	stream() { return stringstream_;}

	std::stringstream stringstream_;
};

class DummyResourceLocator : public ResourceLocator
{
public:
	virtual ~DummyResourceLocator() {}

	/// @brief Create a ResourceStream object from the given resource
	/// source, so that its stream can be passed to the ResourceLoader
	virtual
	ResourceStreamOP
	locate_resource_stream(
		std::string const &
	) const{
		return ResourceStreamOP( new DummyResourceStream );
	}

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP
	) {}

	void
	show(std::ostream & out) const { out << "DummyResourceLocator" << std::endl; }

	std::string
	type() const { return "DummyResourceLocator"; }

};

class DummyResourceLocatorCreator : public ResourceLocatorCreator
{
public:
	virtual
	~DummyResourceLocatorCreator() {}

	virtual
	ResourceLocatorOP
	create_resource_locator() const {
		return ResourceLocatorOP( new DummyResourceLocator );
	}

	virtual
	std::string locator_type() const {
		return "DummyResourceLocator";
	}

};


class DummyResourceFallbackConfigurationBase : public FallbackConfiguration
{
public:

	virtual
	LoaderType
	get_resource_loader( ResourceDescription const & ) const { return "DummyResource"; }

	virtual
	LocatorID
	get_locator_id( ResourceDescription const & ) const { return "DummyResourceLocator"; }

	virtual
	ResourceOptionsOP
	get_resource_options( ResourceDescription const & ) const { return ResourceOptionsOP( new DummyResourceOptions(2) ); }

};

class DummyResourceFallbackConfiguration1 : public DummyResourceFallbackConfigurationBase
{
public:
	virtual bool fallback_specified( ResourceDescription const & ) const { return true; }
	virtual std::string could_not_create_resource_error_message( ResourceDescription const & ) const { return ""; }

};

class DummyResourceFallbackConfiguration2 : public DummyResourceFallbackConfigurationBase
{
public:
	virtual bool fallback_specified( ResourceDescription const & ) const { return false; }
	virtual std::string could_not_create_resource_error_message( ResourceDescription const & ) const { return "test error message"; }
};


class DummyResourceFallbackConfiguration1Creator : public FallbackConfigurationCreator
{
public:
	virtual FallbackConfigurationOP create_fallback_configuration() const { return FallbackConfigurationOP( new DummyResourceFallbackConfiguration1 ); }
	virtual std::string resource_description() const { return "DummyResource1"; }
};

class DummyResourceFallbackConfiguration2Creator : public FallbackConfigurationCreator
{
public:
	virtual FallbackConfigurationOP create_fallback_configuration() const { return FallbackConfigurationOP( new DummyResourceFallbackConfiguration2 ); }
	virtual std::string resource_description() const { return "DummyResource2"; }
};



class JD2ResourceManagerTests : public CxxTest::TestSuite {

private:
	basic::resource_manager::FallbackConfigurationRegistrator< DummyResourceFallbackConfiguration1Creator > reg_fallback_dummy1;
	basic::resource_manager::FallbackConfigurationRegistrator< DummyResourceFallbackConfiguration2Creator > reg_fallback_dummy2;

public:

	void setUp() {
		protocols_init();
	}

	// @brief test default options and default locator
	void test_JD2ResourceManager_without_options() {
		ResourceTag rTag = "the_best_loops_file_ever_made";
		ResourceDescription rDesc = "LoopsFile";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = "protocols/util/demo_loopfile.loops";
		my_config.loader_type = "LoopsFile";
		my_config.resource_options_tag = "";

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
		protocols::loops::LoopsFileDataCOP my_loops = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData const > ( my_resource );

		TS_ASSERT( my_loops ); // make sure we got back the right resource type

		protocols::loops::LoopsFileData lfd = *my_loops;
		TS_ASSERT( lfd.size() == 1 );
		TS_ASSERT( lfd[1].start_res().resindex() == 50 );
		TS_ASSERT( lfd[1].start_res().chain() == 'A' );
		TS_ASSERT( lfd[1].start_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].end_res().resindex() == 60 );
		TS_ASSERT( lfd[1].end_res().chain() == 'A' );
		TS_ASSERT( lfd[1].end_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].cutpoint_res().resindex() == 55 );
		TS_ASSERT( lfd[1].cutpoint_res().chain() == 'A' );
		TS_ASSERT( lfd[1].cutpoint_res().insertion_code() == 'A' );

	}

	// @brief test explicitly passing default options and default locator
	void test_JD2ResourceManager_with_default_options() {
		ResourceTag rTag = "feeling_loopy";
		ResourceOptionsTag resource_options_tag = "use_the_standard_options";
		ResourceDescription rDesc = "LoopsFile";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = "protocols/util/demo_loopfile.loops";
		my_config.loader_type = "LoopsFile";
		my_config.resource_options_tag = resource_options_tag;

		protocols::loops::LoopsFileOptionsOP lfo( new protocols::loops::LoopsFileOptions() );

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		lazy_resource_manager->add_resource_options( resource_options_tag, lfo );
		ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
		protocols::loops::LoopsFileDataCOP my_loops = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData const > ( my_resource );

		TS_ASSERT( my_loops ); // make sure we got back the right resource type

		protocols::loops::LoopsFileData lfd = *my_loops;
		TS_ASSERT( lfd.size() == 1 );
		TS_ASSERT( lfd[1].start_res().resindex() == 50 );
		TS_ASSERT( lfd[1].start_res().chain() == 'A' );
		TS_ASSERT( lfd[1].start_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].end_res().resindex() == 60 );
		TS_ASSERT( lfd[1].end_res().chain() == 'A' );
		TS_ASSERT( lfd[1].end_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].cutpoint_res().resindex() == 55 );
		TS_ASSERT( lfd[1].cutpoint_res().chain() == 'A' );
		TS_ASSERT( lfd[1].cutpoint_res().insertion_code() == 'A' );

	}

	// @brief test explicitly passing nondefault options and default locator
	void test_JD2ResourceManager_with_nondefault_options() {
		ResourceTag rTag = "froot_loops";
		ResourceOptionsTag resource_options_tag = "use_custom_options";
		ResourceDescription rDesc = "LoopsFile";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = "protocols/util/demo_loopfile.loops";
		my_config.loader_type = "LoopsFile";
		my_config.resource_options_tag = resource_options_tag;

		protocols::loops::LoopsFileOptionsOP lfo( new protocols::loops::LoopsFileOptions() );
		lfo->prohibit_single_residue_loops( false );

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		lazy_resource_manager->add_resource_options( resource_options_tag, lfo );
		ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
		protocols::loops::LoopsFileDataCOP my_loops = utility::pointer::dynamic_pointer_cast< protocols::loops::LoopsFileData const > ( my_resource );

		TS_ASSERT( my_loops ); // make sure we got back the right resource type

		protocols::loops::LoopsFileData lfd = *my_loops;
		TS_ASSERT( lfd.size() == 1 );
		TS_ASSERT( lfd[1].start_res().resindex() == 50 );
		TS_ASSERT( lfd[1].start_res().chain() == 'A' );
		TS_ASSERT( lfd[1].start_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].end_res().resindex() == 60 );
		TS_ASSERT( lfd[1].end_res().chain() == 'A' );
		TS_ASSERT( lfd[1].end_res().insertion_code() == 'A' );

		TS_ASSERT( lfd[1].cutpoint_res().resindex() == 55 );
		TS_ASSERT( lfd[1].cutpoint_res().chain() == 'A' );
		TS_ASSERT( lfd[1].cutpoint_res().insertion_code() == 'A' );

	}

	// @brief test omitting necessary data and make sure the JD2ResourceManager throws an exception
	void test_JD2ResourceManager_for_handling_missing_input() {
		ResourceTag rTag = "this_should_not_work";
		ResourceDescription rDesc = "LoopsFile";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = "protocols/util/demo_loopfile.loops";
		my_config.loader_type = "LoopsFile";
		my_config.resource_options_tag = "options_tag_without_a_corresponding_object";

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );

		try {
			ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "Unable to find resource options for the resource options tag 'options_tag_without_a_corresponding_object'.";
			TS_ASSERT( expected_error_message == e.msg() );
		}
	}

	// @brief test making an error and make sure the JD2ResourceManager throws an exception
	void test_JD2ResourceManager_for_handling_wonky_input() {
		ResourceTag rTag = "and_now_for_something_completely_different";
		ResourceOptionsTag resource_options_tag = "use_wrong_type_of_options";
		ResourceDescription rDesc = "LoopsFile";
		JobTag jTag = "fisher_prices_my_first_job";

		ResourceConfiguration my_config;
		my_config.resource_tag = rTag;
		my_config.locator_tag = "";
		my_config.locator_id = "protocols/util/demo_loopfile.loops";
		my_config.loader_type = "LoopsFile";
		my_config.resource_options_tag = resource_options_tag;

		DummyResourceOptionsOP dro( new DummyResourceOptions() );

		ResourceManager * resource_manager = ResourceManager::get_instance();
		LazyResourceManager * lazy_resource_manager = dynamic_cast< LazyResourceManager * > ( resource_manager );

		TS_ASSERT( lazy_resource_manager ); // make sure we got back the right resource_manager type

		lazy_resource_manager->add_resource_tag_by_job_tag( rDesc, jTag, rTag );
		lazy_resource_manager->add_resource_configuration( rTag, my_config );
		lazy_resource_manager->add_resource_options( resource_options_tag, dro );

		try {
			ResourceOP my_resource = lazy_resource_manager->get_resource_by_job_tag( rDesc, jTag );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_error_message = "LoopsFileLoader expected to be given a LoopsFileOptions object, but was given a non-LoopsFileOptions object of type 'DummyResourceOptions', which has the name 'unnamed'.";
			TS_ASSERT( expected_error_message == e.msg() );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_options() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		ResourceOptionsFactory * factory = ResourceOptionsFactory::get_instance();
		factory->factory_register( ResourceOptionsCreatorOP( new DummyResourceOptionsCreator ) );

		std::string xmlfile =
			"<ResourceOptions>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=5/>\n"
			"</ResourceOptions>\n";
		std::istringstream resource_options_stream( xmlfile );
		utility::tag::TagCOP resource_options_tags = utility::tag::Tag::create( resource_options_stream );
		try {
			jd2rm->read_resource_options_tags( resource_options_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		try {
			TS_ASSERT( jd2rm->has_resource_options( "dummyopt1" ) );
			TS_ASSERT( dynamic_cast< DummyResourceOptions * > ( jd2rm->find_resource_options( "dummyopt1" ).get() ));
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	void test_JD2ResourceManager_read_two_dummy_resource_options() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<ResourceOptions>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=5/>\n"
			"  <DummyResourceOptions tag=dummyopt2 somevar=6/>\n"
			"</ResourceOptions>\n";
		std::istringstream resource_options_stream( xmlfile );
		utility::tag::TagCOP resource_options_tags = utility::tag::Tag::create( resource_options_stream );
		try {
			jd2rm->read_resource_options_tags( resource_options_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		try {
			TS_ASSERT( jd2rm->has_resource_options( "dummyopt1" ) );
			TS_ASSERT( dynamic_cast< DummyResourceOptions * > ( jd2rm->find_resource_options( "dummyopt1" ).get() ));
			TS_ASSERT( dynamic_cast< DummyResourceOptions * > ( jd2rm->find_resource_options( "dummyopt1" ).get() )->somevar_ == 5 );
			TS_ASSERT( jd2rm->has_resource_options( "dummyopt2" ) );
			TS_ASSERT( dynamic_cast< DummyResourceOptions * > ( jd2rm->find_resource_options( "dummyopt2" ).get() ));
			TS_ASSERT( dynamic_cast< DummyResourceOptions * > ( jd2rm->find_resource_options( "dummyopt2" ).get() )->somevar_ == 6 );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}


	void test_JD2ResourceManager_read_dummy_resource_options_reused_tag() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<ResourceOptions>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=5/>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=6/>\n"
			"</ResourceOptions>\n";
		std::istringstream resource_options_stream( xmlfile );
		utility::tag::TagCOP resource_options_tags = utility::tag::Tag::create( resource_options_stream );
		try {
			jd2rm->read_resource_options_tags( resource_options_tags );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Duplicated tag, 'dummyopt1' assigned to a ResourceOptions object with type 'DummyResourceOptions'.\n"
				"Prevous ResourceOptions object with this tag was of type 'DummyResourceOptions'\n";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();
		ResourceLoaderFactory::get_instance()->factory_register( ResourceLoaderCreatorOP( new DummyResourceLoaderCreator ) );

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locatorID=1ten.pdb />\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		try {
			TS_ASSERT( jd2rm->has_resource_configuration( "dummyresource" ) );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_missing_locatorID() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource />\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource subtag with LoaderType 'DummyResource' was not supplied with a locatorID tag, which is required.\n";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_duplicate_tag() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locatorID=1ten.pdb/>\n"
			"  <DummyResource tag=dummyresource locatorID=1teb.pdb/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Duplicated tag, 'dummyresource' assigned to a Resource object with ResourceLoader type 'DummyResource'.\n";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_w_resource_options() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string options_xmlfile =
			"<ResourceOptions>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=5/>\n"
			"</ResourceOptions>\n";
		std::istringstream resource_options_stream( options_xmlfile );
		utility::tag::TagCOP resource_options_tags = utility::tag::Tag::create( resource_options_stream );
		jd2rm->read_resource_options_tags( resource_options_tags );

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locatorID=1ten.pdb options=dummyopt1/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cout << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		try {
			TS_ASSERT( jd2rm->has_resource_configuration( "dummyresource" ) );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	void test_JD2ResourceManager_read_dummy_resource_w_resource_options_missing() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locatorID=1ten.pdb options=dummyopt1/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource 'dummyresource' with LoaderType 'DummyResource' was given the tag for a ResourceLoaderOptions, 'dummyopt1', which has not previously been declared.";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_w_file_tag() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource file=1ten.pdb options=dummyopt1/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource 'dummyresource' with LoaderType 'DummyResource' was given the tag for a ResourceLoaderOptions, 'dummyopt1', which has not previously been declared.";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_w_file_and_locatorID_tags() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource file=1ten.pdb locatorID=1ten.pdb options=dummyopt1/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource subtag with LoaderType 'DummyResource' has both file='1ten.pdb' and locatorID='1ten.pdb' but it is only allowed to have one.\n";
			if ( e.msg() != expected_error ) {
				std::cout << "'" << e.msg() << "'" << std::endl;
				std::cout << "'" << expected_error << "'" << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_w_non_file_system_resource_locator_and_file_tags() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locator=\"\" file=1ten.pdb/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource 'dummyresource' with LoaderType 'DummyResource' has locator=\"\" and file='1ten.pdb'. But specifying a file only compatible with the FileSystemResourceLocator.";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_w_file_system_resource_locator_and_file_tags() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locator=FileSystemResourceLocator file=1ten.pdb options=dummyopt1/>\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::string expected_error =
				"Resource subtag with LoaderType 'DummyResource' was given the ResourceLocator tag 'FileSystemResourceLocator', which has not previously been declared.";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	void test_JD2ResourceManager_read_dummy_resource_locator() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();
		ResourceLocatorFactory::get_instance()->factory_register( ResourceLocatorCreatorOP( new DummyResourceLocatorCreator ) );

		std::string xmlfile =
			"<ResourceLocators>\n"
			"  <DummyResourceLocator tag=dummyresloc />\n"
			"</ResourceLocators>\n";
		std::istringstream resource_locators_stream( xmlfile );
		utility::tag::TagCOP resource_locators_tags = utility::tag::Tag::create( resource_locators_stream );
		try {
			jd2rm->read_resource_locators_tags( resource_locators_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( jd2rm->has_resource_locator( "dummyresloc" ) );
	}


	void test_JD2ResourceManager_read_dummy_resource_w_resource_options_and_locator() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		std::string locator_xmlfile =
			"<ResourceLocators>\n"
			"  <DummyResourceLocator tag=dummyresloc />\n"
			"</ResourceLocators>\n";
		std::istringstream resource_locators_stream( locator_xmlfile );
		utility::tag::TagCOP resource_locators_tags = utility::tag::Tag::create( resource_locators_stream );
		jd2rm->read_resource_locators_tags( resource_locators_tags );

		std::string options_xmlfile =
			"<ResourceOptions>\n"
			"  <DummyResourceOptions tag=dummyopt1 somevar=5/>\n"
			"</ResourceOptions>\n";
		std::istringstream resource_options_stream( options_xmlfile );
		utility::tag::TagCOP resource_options_tags = utility::tag::Tag::create( resource_options_stream );
		jd2rm->read_resource_options_tags( resource_options_tags );

		std::string xmlfile =
			"<Resources>\n"
			"  <DummyResource tag=dummyresource locatorID=1ten.pdb options=dummyopt1 locator=dummyresloc />\n"
			"</Resources>\n";
		std::istringstream resources_stream( xmlfile );
		utility::tag::TagCOP resources_tags = utility::tag::Tag::create( resources_stream );
		try {
			jd2rm->read_resources_tags( resources_tags );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			TS_ASSERT( false );
			std::cout << "Raised exception " << e.msg() << std::endl;
		}
		try {
			TS_ASSERT( jd2rm->has_resource_configuration( "dummyresource" ) );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::cerr << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}

	}

	/// Try to instantiate a resource that has not been provided to the JD2ResourceManager
	/// but for which a FallbackConfiguration has been registered with the FallbackConfigurationFactory.
	/// The DummyResource1 resource description can be instantiated, illustrating a case
	/// where the correct set of command line options have been given for a desired resource
	void do_not_test_JD2ResourceManager_fallback_1() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		// This call will fall back on the resource
		TS_ASSERT( jd2rm->has_resource_with_description( "DummyResource1" ));
		ResourceOP dummy;
		try {
			dummy = jd2rm->get_resource( "DummyResource1" );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::cerr  << "Raised exception " << e.msg() << std::endl;
			TS_ASSERT( false );
		}
		TS_ASSERT( dummy );
		DummyResource * dummy_downcasted = dynamic_cast< DummyResource * > ( dummy.get() );
		TS_ASSERT( dummy_downcasted );
		// make sure the ResourceOptions specified by the FallbackConfiguration was used to
		// construct this resource
		TS_ASSERT( dummy_downcasted->somevar_ == 2 );
	}

	/// Try to instantiate a resource that has not been provided to the JD2ResourceManager
	/// but for which a FallbackConfiguration has been registered with the FallbackConfigurationFactory.
	/// The DummyResource2 resource description cannot be instantiated, illustrating a case
	/// where the correct set of command line options have NOT been given for a desired resource
	void test_JD2ResourceManager_fallback_2() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		// This call will fall back on the DummyResourceFallbackConfiguration2
		TS_ASSERT( ! jd2rm->has_resource_with_description( "DummyResource2" ));
		ResourceOP dummy;
		try {
			dummy = jd2rm->get_resource( "DummyResource2" );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_error =
				"JD2ResourceManager does not have a resource corresponding to the resource description 'DummyResource2'. for job 'EMPTY_JOB_use_jd2'.\n"
				"Resources may be specified on the command line or through the JD2ResourceManagerJobInputter resource-declaration file.\n"
				"The FallbackConfiguration for this resource description gives this error:\n"
				"test error message\n"
				"Thrown from JD2ResourceManager::get_resource\n";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}

	/// Try to instantiate a resource that has not been provided to the JD2ResourceManager
	/// and for which no FallbackConfiguration has been registered with the FallbackConfigurationFactory.
	/// This illustrates a case in which the protocol has requested a resource which either
	/// the protocol contains a typo and is requesting one resource when it means another,
	/// or when a resource to be used must be provided through the JD2ResourceManagerJobInputter.
	void test_JD2ResourceManager_fallback_3() {
		JD2ResourceManager * jd2rm = JD2ResourceManager::get_jd2_resource_manager_instance();
		jd2rm->clear();

		// There is no fallback for dummy resource 3
		TS_ASSERT( ! jd2rm->has_resource_with_description( "DummyResource3" ));
		ResourceOP dummy;
		try {
			dummy = jd2rm->get_resource( "DummyResource3" );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_error =
				"JD2ResourceManager does not have a resource corresponding to the resource description 'DummyResource3'. for job 'EMPTY_JOB_use_jd2'.\n"
				"Resources may be specified on the command line or through the JD2ResourceManagerJobInputter resource-declaration file.\n"
				"This resource description does not have a FallbackConfiguration defined.\n"
				"Thrown from JD2ResourceManager::get_resource\n";
			if ( e.msg() != expected_error ) {
				std::cout << e.msg() << std::endl;
				std::cout << expected_error << std::endl;
			}
			TS_ASSERT( e.msg() == expected_error );
		}
	}


};
