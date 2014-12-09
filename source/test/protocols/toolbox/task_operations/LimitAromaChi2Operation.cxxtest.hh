// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/LimitAromaChi2Operation.cxxtest.hh
/// @brief  test suite for LimitAromaChi2Operation
/// @author Nobuasu Koga

// Test headers
#include <cxxtest/TestSuite.h>

#include "platform/types.hh"

#include <test/core/init_util.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/graph/Graph.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <string>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.toolbox.task_operations.LimitAromaChi2Operation.cxxtest");

using namespace core;

class LimitAromaChi2OperationTests : public CxxTest::TestSuite{


public: //typedef


	typedef core::conformation::Residue Residue;
	typedef core::conformation::ResidueCOP ResidueCOP;
	typedef core::graph::GraphOP GraphOP;
	typedef core::pack::task::ResidueLevelTask ResidueLevelTask;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::TaskFactory TaskFactory;
	typedef core::pack::rotamer_set::RotamerSet RotamerSet;
	typedef core::pack::rotamer_set::RotamerSetOP RotamerSetOP;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::toolbox::task_operations::LimitAromaChi2Operation LimitAromaChi2Operation;
	typedef protocols::toolbox::task_operations::LimitAromaChi2OperationOP LimitAromaChi2OperationOP;


public: //data 
	
	
	chemical::ResidueTypeSetCAP residue_set;
	
	
public: //


	LimitAromaChi2OperationTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-ex1aro -ex2aro" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// ------------------------------------------ //	
	RotamerSetOP create_rotset( Size const resid, Pose const & pose, ScoreFunctionOP const scfxn, PackerTaskOP const task ) const
	{
		core::pack::rotamer_set::RotamerSetFactory rsf;
		GraphOP packer_neighbor_graph = core::pack::create_packer_graph( pose, *scfxn, task );
		
		Residue res = pose.residue( resid );
	 	RotamerSetOP rotset = rsf.create_rotamer_set( res );
		rotset->set_resid( resid );
		rotset->build_rotamers( pose, *scfxn, *task, packer_neighbor_graph );
		
		return rotset;
	}
	
	void test_limit_aroma_chi2()
	{
		// read pose
		Pose pose;
		core::import_pose::pose_from_pdb( pose, "protocols/moves/test_in.pdb" );
		pose.update_residue_neighbors();
		
		// the residue number for this test is 10		
		Size resid( 10 ); 
		
		// create scorefxn
		ScoreFunctionOP scfxn = core::scoring::get_score_function();

		// set LimitAromaChi2Operation		
		LimitAromaChi2OperationOP limitar( new LimitAromaChi2Operation );
		limitar->chi2max( 100 );
		limitar->chi2min(  50 );
		limitar->include_trp( false );

		// test default task
		PackerTaskOP dtask( TaskFactory::create_packer_task( pose ) );		
		RotamerSetOP drotset = create_rotset( resid, pose, scfxn, dtask );
		
		Size total( 0 );
		for( Size ii=1; ii<=drotset->num_rotamers(); ii++ )
		{			
			ResidueCOP rop = drotset->rotamer( ii );
			
			if( rop->aa() == core::chemical::aa_tyr ||
				  rop->aa() == core::chemical::aa_phe ||
				  rop->aa() == core::chemical::aa_his ||
				 (rop->aa() == core::chemical::aa_trp && limitar->include_trp() ) )
			{
				utility::vector1< Real > chi( rop->chi() );
				if( chi[ 2 ] < 0 ) {
					chi[ 2 ] += 180.0;
				}
				if( chi[ 2 ] > limitar->chi2max() || chi[ 2 ] < limitar->chi2min() ) {
					total++;					
				}
			}
		}		
		
		// test LimitAromaChi2Operation
		TaskFactory taskf;
		taskf.push_back( limitar );
		PackerTaskOP ptask( taskf.create_task_and_apply_taskoperations( pose ) );	
		RotamerSetOP rotset = create_rotset( resid, pose, scfxn, ptask );
					
		for( Size ii=1; ii<=rotset->num_rotamers(); ii++ )
		{			
			ResidueCOP rop = rotset->rotamer( ii );			
			
			if( rop->aa() == core::chemical::aa_tyr ||
				  rop->aa() == core::chemical::aa_phe ||
				  rop->aa() == core::chemical::aa_his ||
				 (rop->aa() == core::chemical::aa_trp && limitar->include_trp() ) )
			{
				TS_ASSERT( rop->nchi() >= 2 );
				utility::vector1< Real > chi( rop->chi() );
				if( chi[ 2 ] < 0 ) {
					chi[ 2 ] += 180.0;
				}
				
				//	std::cout << ii << " " << rop->name3() << " " << chi[ 1 ] << " " << chi[ 2 ] << " " << chi.size() << std::endl;
				TS_ASSERT( chi[ 2 ] <= limitar->chi2max() && chi[ 2 ] >= limitar->chi2min() );
			}
		}

		// std::cout << drotset->num_rotamers() << " " << rotset->num_rotamers() << std::endl;
		TS_ASSERT( drotset->num_rotamers() == ( rotset->num_rotamers() + total ) );
		
		// test LimitAromaChi2Operation without TRP
		taskf.clear();
		limitar->include_trp( true );
		taskf.push_back( limitar );
		ptask = taskf.create_task_and_apply_taskoperations( pose );	
		rotset = create_rotset( resid, pose, scfxn, ptask );
		
		for( Size ii=1; ii<=rotset->num_rotamers(); ii++ )
		{			
			ResidueCOP rop = rotset->rotamer( ii );
	
			if( rop->aa() == core::chemical::aa_tyr ||
					rop->aa() == core::chemical::aa_phe ||
					rop->aa() == core::chemical::aa_his ||
				 (rop->aa() == core::chemical::aa_trp && limitar->include_trp() ) )
				{
				TS_ASSERT( rop->nchi() >= 2 );
				utility::vector1< Real > chi( rop->chi() );
				if( chi[ 2 ] < 0 ) {
					chi[ 2 ] += 180.0;
				}
				
				// std::cout << ii << " " << rop->name3() << " " << chi[ 1 ] << " " << chi[ 2 ] << " " << chi.size() << std::endl;
				TS_ASSERT( chi[ 2 ] <= limitar->chi2max() && chi[ 2 ] >= limitar->chi2min() );
			}
		}
	}


};

