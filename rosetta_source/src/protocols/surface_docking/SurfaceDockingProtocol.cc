// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
//     vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file SurfaceDockingProtocol.cc
/// @author Robin A Thottungal (raugust1@jhu.edu)

// Unit Headers
#include <protocols/surface_docking/SurfaceDockingProtocol.hh>

// Package Headers
#include <protocols/surface_docking/CentroidRelaxMover.hh>
#include <protocols/surface_docking/FullatomRelaxMover.hh>
#include <protocols/surface_docking/SurfaceOrientMover.hh>
//#include <protocols/surfaceDocking/SurfaceOrientMover.fwd.hh>

// Project headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/moves/MinMover.fwd.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>

//Utility Headers
#include <utility/exit.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <core/kinematics/FoldTree.hh>
#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;

using namespace protocols::surface_docking;
using namespace protocols;
using namespace protocols::moves;
//using namespace ObjexxFCL::fmt;

////using core::pose::datacache::CacheableDataType::SURFACE_PARAMS;
static basic::Tracer TR("protocols.SurfaceDocking.SurfaceDockingProtocol");

namespace protocols {
namespace surface_docking {

using namespace core;
using protocols::jd2::JobDistributor;

//constructor
SurfaceDockingProtocol::SurfaceDockingProtocol() : Mover(){
	Mover::type( "SurfaceDockingProtocol");
	}

//destructor
SurfaceDockingProtocol::~SurfaceDockingProtocol() {}

void SurfaceDockingProtocol::setupFoldTree(pose::Pose & pose){
	// Make sure the HETATM is first
	TR<<"Total Residue:"<<pose.total_residue()<<std::endl;
	kinematics::FoldTree ft=pose.fold_tree();
	//TR<<"Old FoldTree"<<ft<<std::endl;
	TR<<"FoldTree"<<ft<<std::endl;
	TR<<"Number of Jumps:"<<pose.num_jump()<<std::endl;
}


void SurfaceDockingProtocol::apply(pose::Pose & pose){

	//AddPyMolObserver(pose,false);
	PyMolMoverOP my_PymolMover=new PyMolMover();
	pose.dump_pdb("PymolObserver_testing.pdb");
	my_PymolMover->apply(pose);
	// Setting up the foldTree
	setupFoldTree(pose);

	simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	simple_moves::SwitchResidueTypeSetMover to_all_atom( "fa_standard" );

	//Applying surfaceOrient Mover
	surface_docking::SurfaceOrientMoverOP sf=
                                new surface_docking::SurfaceOrientMover();
	sf->apply(pose);

	//split the pose into separate chains
	utility::vector1< pose::PoseOP > singlechain_poses;
	singlechain_poses = pose.split_by_chain();

	// singlechain_poses[2] contains the protein; can be made generic
	// for multiple proteins too
	to_centroid.apply( *singlechain_poses[2] );
	surface_docking::CentroidRelaxMoverOP centroidrelax =
	                            new surface_docking::CentroidRelaxMover();
	centroidrelax->set_nmoves(10);
	centroidrelax->apply(*singlechain_poses[2]);
	to_all_atom.apply(*singlechain_poses[2]);
	/**
	// joining the pose together
	Size cutpoint ( (*singlechain_poses[1]).total_residue() );
	for ( Size i=1; i<=(*singlechain_poses[2]).total_residue(); ++i ) {
		conformation::ResidueCOP new_rsd = (*singlechain_poses[2]).residue(i).clone();
		if ( i == 1 ) {
			// First residue, connecting with the previous segment through
			// a jump
			(*singlechain_poses[1]).append_residue_by_jump( *new_rsd, cutpoint );
		}
		else {
			(*singlechain_poses[1]).append_residue_by_bond( *new_rsd );
		}
	}
	pose=*singlechain_poses[1];
	**/
	//pose.copy_segment( tmp_pose.total_residue(), tmp_pose, 1, pose.total_residue()+1 );
	//core::pose::PoseOP new_pose;
	pose.copy_segment((*singlechain_poses[2]).total_residue(),
			 (*singlechain_poses[2]),
			 (*singlechain_poses[1]).total_residue()+1,1);
	setupFoldTree(pose);
	//AddPyMolObserver(pose,false);
	//Fullatom Relax starts!

	//to_all_atom.apply( pose );
	surface_docking::FullatomRelaxMoverOP allatomrelax =
                                new surface_docking::FullatomRelaxMover();
	//allatomrelax->set_nmoves(10);
	core::Size lj_ramp_cycle=5; // inside the constructor
	// set the random number
	core::Real lj_increment = ( 1.0 - 0.02 )/ lj_ramp_cycle;
	for (Size i=1;i<=lj_ramp_cycle;++i){
		allatomrelax->set_smallmovesize(30/i);
		allatomrelax->set_ljrepulsion_weight(0.02+i*lj_increment);
		allatomrelax->set_ecounter(i); // when value of i matches the random number
                      // in the allatomrelax, protein is slide into the surface
		allatomrelax->apply(pose);
	}
	// Final Side Chain re-packing using rtmin; adopted from Dave's protocol


	//core::scoring::SurfaceParameters & surfaceVectors=
	//	*( static_cast< core::scoring::SurfaceParameters * >
	//				( pose.data().get_ptr( SURFACE_PARAMS )() ));

	SurfaceParametersOP surfaceVectors= new SurfaceParameters();

	// Creating a job compatible with JD2
	static protocols::jd2::JobOP job
                        =jd2::JobDistributor::get_instance()->current_job();
	std::string job_name (JobDistributor::get_instance()->
                                       job_outputter()->output_name( job ) );

	job->add_string_string_pair(surfaceVectors->strSURFA0," ");
	job->add_string_string_pair(surfaceVectors->strSURFA1," ");
	job->add_string_string_pair(surfaceVectors->strSURFA2," ");
	JobDistributor::get_instance()->job_outputter()->
                                          other_pose( job,pose, "Surface_");
	JobDistributor::get_instance()->job_outputter()->
						other_pose( job,*singlechain_poses[1], "Surface-1_");
	JobDistributor::get_instance()->job_outputter()->
                       	other_pose( job,*singlechain_poses[2], "Surface-2_");

}

std::string SurfaceDockingProtocol::get_name() const {
	return "SurfaceDockingProtocol";
	}

}	//surfaceDockingProtocol

}	//protocol
