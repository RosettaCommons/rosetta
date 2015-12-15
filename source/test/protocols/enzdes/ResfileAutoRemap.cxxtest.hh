// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/ResfileAutoRemap.cxxtest.hh
/// @brief  test suite for making resfiles work with poses changing length
/// @author Florian Richter

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/protocols/init_util.hh>

#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh> //need for additional residue
#include <core/chemical/ResidueTypeSet.hh>

#include <basic/options/option.hh> //needed to set option
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //function for reading cstfiles
#include <protocols/enzdes/EnzdesTaskOperations.hh>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.enzdes.ResfileAutoRemap.cxxtest");

using namespace core;

/// @detail test that checks whether a resfile can be applied to a
/// pose that had its length changed.
class ResfileAutoRemapTest : public CxxTest::TestSuite
{

public:
	ResfileAutoRemapTest() {};
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP enz_io;


	// Shared initialization goes here.
	void setUp() {
		protocols_init();
		// Residue definitions can't be supplied on the command line b/c
		// the ResidueTypeSet is already initialized.
		using namespace core::chemical;
		utility::vector1< std::string > params_files;
		ResidueTypeSetCOP const_residue_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		ResidueTypeSet & residue_set = const_cast< ResidueTypeSet & >(*const_residue_set);
		if ( !residue_set.has_name("D2N") ) params_files.push_back("protocols/enzdes/D2N.params");
		residue_set.read_files_for_custom_residue_types(params_files);
		basic::options::option[basic::options::OptionKeys::run::preserve_header ].value(true);

		enz_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO(const_residue_set) );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_resfile_auto_remap()
	{
		using namespace core::pose;
		using namespace core::pose::datacache;
		using namespace core::pack::task::operation;
		using namespace core::scoring::constraints;

		//typedef core::id::AtomID AtomID;
		bool optionsaveval =  basic::options::option[basic::options::OptionKeys::enzdes::detect_design_interface ].value();
		basic::options::option[basic::options::OptionKeys::enzdes::detect_design_interface ].value(true);

		pose::Pose test_pose;
		core::pack::task::TaskFactory task_factory, task_factory_compare;

		core::import_pose::pose_from_pdb( test_pose, "protocols/enzdes/ligtest_it.pdb");
		scoring::ScoreFunctionOP scorefxn = scoring::ScoreFunctionFactory::create_score_function("enzdes");

		//now let's use the enzdes machinery to read in a cstfile and generate
		//the constraint set, results should be identical to manually created constraints
		enz_io->read_enzyme_cstfile("protocols/enzdes/ligtest_it.cst");
		enz_io->add_constraints_to_pose(test_pose, scorefxn, false);

		core::pose::Pose compare_pose = test_pose;

		(*scorefxn)(compare_pose);
		(*scorefxn)(test_pose);

		//muck around with the pose
		test_pose.observer_cache().set( core::pose::datacache::CacheableObserverType::LENGTH_EVENT_COLLECTOR, CacheableObserverOP( new core::pose::datacache::LengthEventCollector() ) );
		test_pose.delete_polymer_residue( 67 );
		test_pose.delete_polymer_residue( 78 );
		(*scorefxn)(test_pose);
		//set up the task operations
		//protocols::enzdes::DetectProteinLigandInterfaceOP interface_op = new protocols::enzdes::DetectProteinLigandInterface;

		//task_factory.push_back( interface_op );
		task_factory.push_back( TaskOperationCOP( new protocols::enzdes::SetCatalyticResPackBehavior ) );
		task_factory.push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfileAndObeyLengthEvents( "protocols/enzdes/resfile_remap.resfile") ) );

		//task_factory_compare.push_back( interface_op );
		task_factory_compare.push_back( TaskOperationCOP( new protocols::enzdes::SetCatalyticResPackBehavior ) );
		task_factory_compare.push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile( "protocols/enzdes/resfile_remap.resfile") ) );
		//have to wipe out pdbinfo, this test is testing non-pdbinfo remapping functionality
		test_pose.pdb_info( PDBInfoOP( new core::pose::PDBInfo( test_pose ) ));
		core::pack::task::PackerTaskOP ptask( task_factory.create_task_and_apply_taskoperations( test_pose ) );
		core::pack::task::PackerTaskOP ptask_compare( task_factory_compare.create_task_and_apply_taskoperations( compare_pose ) );

		//compare tasks a bit
		TS_ASSERT( ptask->residue_task( 97 ).include_current() == ptask_compare->residue_task( 99 ).include_current() );
		TS_ASSERT( ptask->residue_task( 99 ).include_current() != ptask_compare->residue_task( 99 ).include_current() );
		TS_ASSERT( ptask_compare->residue_task( 74 ).being_designed() );
		TS_ASSERT( !ptask->residue_task( 74 ).being_designed() );
		TS_ASSERT( !ptask_compare->residue_task( 73 ).being_designed() );
		TS_ASSERT( ptask->residue_task( 73 ).being_designed() );

		TR << "pos 97 include_current in ptask " << ptask->residue_task( 97 ).include_current() << ", pos 99 " << ptask->residue_task( 99 ).include_current() << std::endl;

		//muck around with the pose a little more
		core::conformation::ResidueOP dummyres = test_pose.residue(3).clone();
		test_pose.append_polymer_residue_after_seqpos( *dummyres, 19, true );
		test_pose.append_polymer_residue_after_seqpos( *dummyres, 57, true );
		test_pose.append_polymer_residue_after_seqpos( *dummyres, 95, true );
		(*scorefxn)(test_pose);
		test_pose.pdb_info( PDBInfoOP( new core::pose::PDBInfo( test_pose ) ));

		//compare the tasks a bit more
		core::pack::task::PackerTaskOP ptask2( task_factory.create_task_and_apply_taskoperations( test_pose ) );

		TS_ASSERT( ptask2->residue_task( 99 ).include_current() != ptask_compare->residue_task( 99 ).include_current() );
		TS_ASSERT( ptask2->residue_task( 100 ).include_current() == ptask_compare->residue_task( 99 ).include_current() );
		TS_ASSERT( !ptask2->residue_task( 74 ).being_designed() );
		TS_ASSERT( ptask2->residue_task( 75 ).being_designed() );


		TR << "pos 97 include_current in ptask " << ptask2->residue_task( 97 ).include_current() << ", pos 99 " << ptask2->residue_task( 99 ).include_current() << ", pos 100 " << ptask2->residue_task( 100 ).include_current() << std::endl;


		basic::options::option[basic::options::OptionKeys::enzdes::detect_design_interface ].value(optionsaveval);
	}


};

