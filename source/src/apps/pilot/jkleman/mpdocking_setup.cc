// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mpdocking.cc
/// @brief   Dock two membrane proteins in the membrane
/// @details last Modified: 4/4/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/util.hh>
#include <core/membrane/geometry/util.hh>
#include <core/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
//#include <protocols/membrane/OptimizeMembraneMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace core::membrane::geometry;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::membrane;
using namespace core::membrane::geometry;

static thread_local basic::Tracer TR( "apps.pilot.jkleman.mpdocking_setup" );

////////////////////////////////////////////////////////////////////////////////

// vector show function
template< typename T_ >
void show( utility::vector1< T_ > vector){
	for ( Size i = 1; i <= vector.size(); ++i ){
		TR << utility::to_string(vector[i]) << ", ";
	}
	TR << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

// read poses
utility::vector1< PoseOP > read_poses() {
	
	// cry if PDBs not given
	if ( ! option[OptionKeys::in::file::s].user() ){
		utility_exit_with_message("Please provide two PDB files!");
	}

	// get filenames from Optionsystem
	utility::vector1< std::string > pdbfiles = option[OptionKeys::in::file::s]();

	show( pdbfiles );
	utility::vector1< PoseOP > poses;

	// put poses vector
	for ( Size i = 1; i <= pdbfiles.size(); ++i ){
		
		// read poses
		PoseOP pose = core::import_pose::pose_from_pdb( pdbfiles[i] );
				
		// put poses and topologies into vector
		poses.push_back( pose );
	}
	
	return poses;

}// read poses

////////////////////////////////////////////////////////////////////////////////

// get topologies, poses only needed for nres check
utility::vector1< SpanningTopologyOP > get_topologies( utility::vector1< std::string > spanfiles, utility::vector1< PoseOP > poses ){

	// check for file correspondence
	if ( spanfiles.size() != poses.size() ){
		utility_exit_with_message("Number of PDB files does not match number of spanfiles! They must match!");
	}
		
	// initialize variables
	show( spanfiles );
	utility::vector1< SpanningTopologyOP > topologies;
	
	// put spanfiles into vector
	for ( Size i = 1; i <= spanfiles.size(); ++i ){
		
		TR << "spanfiles" << spanfiles[i] << std::endl;
		TR << "poses" << poses[i]->total_residue() << std::endl;
		
		// read spanfiles
		SpanningTopologyOP topo = new SpanningTopology( spanfiles[i], poses[i]->total_residue() );
		
		// put topologies into vector
		topologies.push_back( topo );
	}

	for ( Size i = 1; i <= topologies.size(); ++i ){
		topologies[i]->show();
	}
		
	// make sure lengths of poses and topologies are the same
	for ( Size i = 1; i <= poses.size(); ++i ){
			
		if ( poses[i]->total_residue() != topologies[i]->nres_topo() ){
			utility_exit_with_message("Number of residues in PDB " + utility::to_string(i) + " doesn't match the number of residues in the corresponding spanfile!");
		}
	}
	
	return topologies;
}// read spanfiles

////////////////////////////////////////////////////////////////////////////////

// script that creates input for MPDocking protocol
// = reads in two PDB files and 2 spanfiles
// = attaches the membrane virtual residue to pose 1
// = sets initial position of pose 1
// = optimizes the membrane position of pose 2
// = rotates and translates pose 2 such that the membranes superimpose
// = concatenates the PDB files
// = concatenates the spanfiles
// = slides together the 2 poses into contact
// = dumps PDB and spanfiles
// => CURRENTLY ONLY WORKS FOR 2 BODY DOCKING! (the slide together move would need to be adapted for multiple docking partners)

void mpdocking_setup(){

	TR << "mpdocking_setup" << std::endl;
	TR << "PLEASE MAKE SURE YOUR PDBS AND SPANFILES CORRESPOND TO EACH OTHER IN THE INPUT VECTORS!!!" << std::endl;

	// read in poses
	utility::vector1< PoseOP > poses( read_poses() );
	
	// read spanfiles
	utility::vector1< std::string > spanfiles( spanfile_names() );

	// get topologies
	utility::vector1< SpanningTopologyOP > topologies( get_topologies( spanfiles, poses ) );

	TR << "read topologies" << std::endl;

	// get poses and topologies from vectors
	PoseOP pose1 = poses[1];
	PoseOP pose2 = poses[2];
	std::string spanfile1 = spanfiles[1];
	std::string spanfile2 = spanfiles[2];
	SpanningTopologyOP topo1 = topologies[1];
	SpanningTopologyOP topo2 = topologies[2];

	TR << "got poses and topologies" << std::endl;

	// OPTIMIZE MEMBRANE OF POSE1
	// 1) default EmbeddingConfig
	EmbeddingDefOP membrane1 = new EmbeddingDef();

	TR << "set EmbedDef" << std::endl;

//	TR << "Pre-add" << std::endl;
//	TR << "res 1: " << pose1->residue( 1 ).name3() << std::endl;
//	pose1->residue(1).atom(2).xyz().show();
//	TR << std::endl;
//	TR << "res 15: " << pose1->residue( 15 ).name3() << std::endl;
//	pose1->residue(15).atom(2).xyz().show();
//	TR << std::endl;
	
	// 2) attach MEM to pose 1
	AddMembraneMoverOP add_membrane1 = new AddMembraneMover( membrane1->center(), membrane1->normal(), spanfile1, true);
	add_membrane1->apply( *pose1 );

//	TR << "post-add" << std::endl;
//	TR << "res 1: " << pose1->residue( 1 ).name3() << std::endl;
//	pose1->residue(1).atom(2).xyz().show();
//	TR << std::endl;
//	TR << "res 15: " << pose1->residue( 15 ).name3() << std::endl;
//	pose1->residue(15).atom(2).xyz().show();
//	TR << std::endl;

	TR << "attached MEM to pose1" << std::endl;
	
//	TR << "before initial" << std::endl;
//	pose1->conformation().show_membrane();
	
	// show foldtree
	pose1->fold_tree().show(std::cout);
	
	// reorder only reorders, but does not rename jump edges
	core::kinematics::FoldTree foldtree = pose1->fold_tree();
	foldtree.reorder( pose1->conformation().membrane_info()->membrane_rsd_num() );
	pose1->fold_tree( foldtree );
	
	// show foldtree
	TR << "foldtree reordered" << std::endl;
	pose1->fold_tree().show(std::cout);

	// before move
	pose1->dump_pdb("before.pdb");

	// 3) set initial position of pose 1 in the membrane
	MembranePositionFromTopologyMoverOP initial_position1 = new MembranePositionFromTopologyMover( true );
	TR << "constructor called" << std::endl;
	initial_position1->apply( *pose1 );

//	TR << "after initial" << std::endl;
//	pose1->conformation().show_membrane();

//	TR << "set initial position of pose1" << std::endl;
	
	// 4) optimize membrane position of pose 1
//	OptimizeMembraneMoverOP optimize_membrane1 = new OptimizeMembraneMover( true, true );
//	optimize_membrane1->apply( *pose1 );

	// OPTIMIZE MEMBRANE OF POSE2
	// 1) copy pose 2 and do the same procedure (add membrane, get initial position
	//    and optimize membrane)
	//	  don't want to do this on pose2 itself because then I can't remove the membrane
	//    VRT any more before concatenating it to pose1
//	PoseOP pose2a( pose2->clone() );
//
//	TR << "copied pose2" << std::endl;
//
//	// 2) add MEM to pose 2a
//	AddMembraneMoverOP add_membrane2a = new AddMembraneMover( membrane1->center(), membrane1->normal(), spanfile2);
//	add_membrane2a->apply( *pose2a );
//
//	TR << "added MEM to pose2" << std::endl;
//	
//	// 3) set initial position of pose 2a in membrane
////	MembranePositionFromTopologyMoverOP initial_position2a = new MembranePositionFromTopologyMover( true );
////	initial_position2a->apply( *pose2a );
//
//	TR << "set initial position of pose2a" << std::endl;
//
//	// 4) optimize membrane position of pose 2a
////	OptimizeMembraneMoverOP optimize_membrane2a = new OptimizeMembraneMover( true, true );
////	optimize_membrane2a->apply( *pose2a );
//	
//	// rotate and translate pose 2a such that membranes superimpose
//	Vector normal1 = pose1->conformation().membrane_center();
//	Vector center1 = pose1->conformation().membrane_normal();
//
//	TR << "got normal and center of pose1" << std::endl;
//
//	SetMembranePositionMoverOP rot_trans = new SetMembranePositionMover( center1, normal1 );
//	rot_trans->apply( *pose2a );
//
//	TR << "rotated and translated" << std::endl;
//	
//	// get membrane position of pose2a
//	Vector normal2a = pose2a->conformation().membrane_center();
//	Vector center2a = pose2a->conformation().membrane_normal();
//	
//	TR << "got normal and center of pose2a" << std::endl;
//
//	// superimpose pose2 with pose2a
//	protocols::toolbox::CA_superimpose( *pose2a, *pose2 );
//
//	TR << "superimposed" << std::endl;
//
//	// make a copy of pose1 to concatenate the other pose to
//	PoseOP pose1a = pose1;
//
//	TR << "copied pose1" << std::endl;
//
//	// concatenate poses (new chain = true)
//	append_pose_to_pose( *pose1a, *pose2, true );
//
//	TR << "appended to poses" << std::endl;
//
//	// concatenate topologies:
//	// 1) copy topo1
//	SpanningTopologyOP topo1a = topo1;
//	
//	// 2) nres in pose1 is total minus membrane residue
//	Size nres_pose1a = pose1a->total_residue() - 1;
//
//	TR << "nres" << std::endl;
//	
//	// 3) add shifted spans from pose2 to topo1a
//	for ( Size i = 1; i <= topo2->nspans(); ++i ){
//		SpanOP shifted_span = topo2->span(i);
//		shifted_span->shift( nres_pose1a );
//		topo1a->add_span( shifted_span );
//	}
//
//	TR << "shifted spans" << std::endl;
//
//	// update MembraneInfo object
//	Size membrane_resnum = pose1a->conformation().membrane_info()->membrane_rsd_num();
//	pose1a->conformation().setup_membrane( membrane_resnum, topo1a );
//
//	TR << "updated MemInfo" << std::endl;
//	
//	// show foldtree
//	pose1a->fold_tree().show(std::cout);
	
//	pull proteins apart by 100 Angstrom
//	Vector trans (100, 0, 0);
//	Vector new_center = center1 + trans;
//	SetMembraneCenterMoverOP translate = new SetMembraneCenterMover( new_center );
//	translate->apply( *pose2 );
	
	// slide together poses into contact:
	// get attractive score
		
	pose1->dump_pdb("pose1after.pdb");	// original pose1
	pose2->dump_pdb("pose2after.pdb");	// original pose2
//	pose1a->dump_pdb("pose1a.pdb"); // full pose
//	pose2a->dump_pdb("pose2a.pdb"); // translated pose2
	topo1->show();				// original topo1
	topo2->show();				// original topo2
//	topo1a->show();				// full topo
	
	TR << "" << std::endl;
	
	
	


}// mpdocking_setup



////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		// initialize option system, random number generators, and all factory-registrators
		devel::init(argc, argv);

		// call my function
		mpdocking_setup();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
