// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rosetta_scripts/RosettaScriptsParser.cxxtest.hh
/// @brief Tests parsing of top-level RosettaScripts blocks that are controlled by DataLoaders.
/// @details Tests that each DataLoader class recognizes its own block and correctly adds its members to the DataMap.
/// Whenever new DataLoader classes are added, they should be tested here.
/// @author Sharon Guffy

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/util/pose_funcs.hh>
// Package headers
#include <protocols/parser/DataLoaderFactory.hh>
#include <protocols/parser/ScoreFunctionLoader.hh>
#include <protocols/parser/FragSetLoader.hh>
#include <protocols/parser/TaskOperationLoader.hh>
#include <protocols/parser/MoveMapFactoryLoader.hh>
#include <protocols/parser/JumpSelectorLoader.hh>
#include <protocols/qsar/scoring_grid/ScoringGridLoader.hh>
#include <protocols/parser/MonteCarloLoader.hh>
#include <protocols/parser/ResidueSelectorLoader.hh>
#include <protocols/parser/ConstraintGeneratorLoader.hh>
#include <protocols/ligand_docking/LigandDockingLoaders.hh>
#include <protocols/loops/loops_definers/LoopsDefinerLoader.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// Utility headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Numeric headers

// C++ headers
#include <string>
#include <iostream>

static THREAD_LOCAL basic::Tracer TR("protocols.jd2.parser.DataLoaders.cxxtest");
////////////////////////////////////////////////////////////////////////
// Tests

class DataLoaderTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP pose_;
public:
	void setUp() {
		protocols_init();
		pose_ = fullatom_poseop_from_string( test_in_pdb_string() );
	}
	void tearDown() {

	}
	//Data loaders:
	///@brief Test that the DataLoaderFactory
	void test_DataLoaderFactory(){
		using namespace protocols::parser;
		///@details The RosettaScripts parser instantiates data loaders by passing the tag names to the DataLoaderFactory
		///This is basically just making sure that all of them are registered under the correct name.
		DataLoaderOP loader;
		//ConstraintGenerator
		std::string name;
		name = ConstraintGeneratorLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//FragSet--looks like this one isn't actually registered
		/*
		name = FragSetLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		*/
		//JumpSelector
		name = JumpSelectorLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//MonteCarlo
		name = MonteCarloLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//MoveMapFactory
		name = MoveMapFactoryLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//ResidueSelector
		name = ResidueSelectorLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//ScoreFunction
		name = ScoreFunctionLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//TaskOperation
		name = TaskOperationLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//ScoringGrid
		name = protocols::qsar::scoring_grid::ScoringGridLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//InterfaceBuilder
		name = protocols::ligand_docking::InterfaceBuilderLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//MoveMapBuilder
		name = protocols::ligand_docking::MoveMapBuilderLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//LigandArea
		name = protocols::ligand_docking::LigandAreaLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
		//LoopsDefiner
		name = protocols::loops::loops_definers::LoopsDefinerLoader::loader_name();
		TS_ASSERT_THROWS_NOTHING( loader = DataLoaderFactory::get_instance()->newDataLoader( name ) );
	}

	//Any simple subclass of these can be used to test the loader.
	//ConstraintGenerator
	void test_ConstraintGeneratorLoader(){
		using namespace protocols::parser;
		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<CONSTRAINT_GENERATORS>\n<ResidueTypeConstraintGenerator name=\"dummy\" />\n</CONSTRAINT_GENERATORS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::ConstraintGeneratorLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "ConstraintGenerators" ) );
		TS_ASSERT( data.has( "ConstraintGenerators", "dummy" ) );
	}
	//ResidueSelector
	void test_ResidueSelectorLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<RESIDUE_SELECTORS>\n<True name=\"dummy\" />\n</RESIDUE_SELECTORS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::ResidueSelectorLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "ResidueSelector" ) );
		TS_ASSERT( data.has( "ResidueSelector", "dummy" ) );
	}
	//ScoreFunction
	void test_ScoreFunctionLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<SCOREFXNS>\n<ScoreFunction name=\"dummy\" />\n</SCOREFXNS>" << std::endl; //use default weights
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::ScoreFunctionLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "scorefxns" ) );
		TS_ASSERT( data.has( "scorefxns", "dummy" ) );
	}
	//TaskOperation
	void test_TaskOperationLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<TASKOPERATIONS>\n<RestrictToRepacking name=\"dummy\" />\n</TASKOPERATIONS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::TaskOperationLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "task_operations" ) );
		TS_ASSERT( data.has( "task_operations", "dummy" ) );
	}
	//JumpSelector
	void test_JumpSelectorLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<JUMP_SELECTORS>\n<Interchain name=\"dummy\" />\n</JUMP_SELECTORS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::JumpSelectorLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "JumpSelector" ) );
		TS_ASSERT( data.has( "JumpSelector", "dummy" ) );
	}

	//MonteCarlo
	void test_MonteCarloLoader(){
		using namespace protocols::parser;

		//All subtags are named MonteCarlo
		basic::datacache::DataMap data;
		std::stringstream ss;
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		data.add( "scorefxns", "dummy2", scorefxn );
		ss << "<MONTECARLOS>\n<MonteCarlo name=\"dummy\" scorefxn=\"dummy2\" />\n</MONTECARLOS>" << std::endl; //Will use default score function
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::MonteCarloLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "montecarlos" ) );
		TS_ASSERT( data.has( "montecarlos", "dummy" ) );
	}
	//MoveMapFactory
	void test_MoveMapFactoryLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<MOVE_MAP_FACTORIES>\n<MoveMapFactory name=\"dummy\" />\n</MOVE_MAP_FACTORIES>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::parser::MoveMapFactoryLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "MoveMapFactory" ) );
		TS_ASSERT( data.has( "MoveMapFactory", "dummy" ) );
	}

	//MoveMapBuilder
	void test_MoveMapBuilderLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<MOVEMAP_BUILDERS>\n<MoveMapBuilder name=\"dummy\" />\n</MOVEMAP_BUILDERS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::ligand_docking::MoveMapBuilderLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "movemap_builders" ) );
		TS_ASSERT( data.has( "movemap_builders", "dummy" ) );
	}
	//LoopsDefiner
	void test_LoopsDefinerLoader(){
		using namespace protocols::parser;

		basic::datacache::DataMap data;
		std::stringstream ss;
		ss << "<LOOP_DEFINITIONS>\n<Loops name=\"dummy\">\n<loop start=\"1\" stop=\"3\" />\n</Loops>\n</LOOP_DEFINITIONS>" << std::endl;
		utility::tag::TagOP tag( new utility::tag::Tag() );
		tag->read( ss );

		protocols::loops::loops_definers::LoopsDefinerLoader loader;
		//Check that the loader can load the tag
		TS_ASSERT_THROWS_NOTHING( loader.load_data( *pose_, tag, data ) );
		//Check that the tag's contents are found in the DataMap
		TS_ASSERT( data.has_type( "loops_definers" ) );
		TS_ASSERT( data.has( "loops_definers", "dummy" ) );
	}

	//TODO: Tests for all the other data loaders
	/*
	//FragSet //complex setup
	void TODOtest_FragSetLoader(){

	//ScoringGrid //This one's a little weird--the outer tag can take options, and frankly I don't understand enough about what a ScoringGrid is or what it's doing to deal with all that (guffysl)
	void TODOtest_ScoringGridLoader(){

	}
	//InterfaceBuilder--these depend on ligand areas, may be hard to test independently
	void TODOtest_InterfaceBuilderLoader(){

	}
	//LigandArea //complex setup
	void TODOtest_LigandAreaLoader(){

	}
	*/
};
