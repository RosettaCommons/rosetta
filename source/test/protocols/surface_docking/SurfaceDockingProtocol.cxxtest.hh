// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/surface_docking/SurfaceVectorLoader.cxxtest.hh
/// @brief test suite for protocols/surface_docking/SurfaceVectorLoader.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/surface_docking/SurfaceVectorLoader.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <protocols/loops/LoopsFileOptions.hh> //this is crazy
#include <protocols/surface_docking/SurfaceDockingProtocol.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/solid_surface/SurfaceEnergies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/surface_docking/SurfaceOrientMover.hh>
#include <basic/Tracer.hh>

// Utility headers


// Numeric headers

// C++ headers
#include <string>

using namespace protocols::surface_docking;

static basic::Tracer TR("surface_orient_mover_test");

class SurfaceDockingProtocolTest : public CxxTest::TestSuite {

public:
	void setUp() {
		protocols_init();
	}

	// @brief test slide_away_from_surface
	void test_slide_away_from_surface() {
		SurfaceVectorLoader loader;
		std::string surf_vec_file( "18.671 -38.353   9.571\n17.705 -36.371   5.096\n12.153 -43.099   8.876" );
		std::string surf_vec_file2( "-28.188  33.809   3.201\n-46.450  30.380   6.325\n-31.129  24.953   4.699\n" );
		std::istringstream lstream( surf_vec_file );
		
		protocols::loops::LoopsFileOptions opts;
		 
		utility::pointer::ReferenceCountOP resource = loader.create_resource( opts, "unit_test", lstream );

		SurfaceParametersOP spptr = utility::pointer::dynamic_pointer_cast< protocols::surface_docking::SurfaceParameters > ( resource );
		
		core::pose::Pose pose;
		
		core::import_pose::pose_from_file(pose, "/Users/mpacella/Rosetta_Surface_Test/surf_orient_test_before.pdb", core::import_pose::PDB_file);
		
		core::Vector protein_centroid, surf_centroid;
		core::Size const rb_jump=pose.num_jump();
		// Last jump is the one that connect the surface with protein
				protocols::geometry::centroids_by_jump (pose, rb_jump,surf_centroid,protein_centroid);
		
	
		//spptr->generate_surface_parameters(SurfaceCG, ProteinCG);
		//TR<<"AB surf vec: "<<spptr->vecAB()<<std::endl;
		core::scoring::solid_surface::SurfaceEnergiesOP surfEs( new core::scoring::solid_surface::SurfaceEnergies );
		surfEs->set_total_residue( pose.total_residue() );
		surfEs->set_residue_range_not_surface( pose.num_jump()+1, pose.total_residue() );
		pose.set_new_energies_object( surfEs );
		
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		
		core::Real score_before = (*sfxn)(pose);
		
		TR<<"score_before: "<<score_before<<std::endl;
		
		core::Vector coords_before = pose.residue(pose.total_residue()).xyz("Ca2p");
	
		
		core::Vector const slide_axis = spptr->slide_axis();
		core::Vector slide_into, slide_away;
		
		
		//asses which side of the surface the protein is and pick the correct direction for the sliding vectors accordingly
		if (surf_centroid.distance_squared(protein_centroid + slide_axis) < surf_centroid.distance_squared(protein_centroid))
		{
			slide_into = spptr->slide_axis();
		}
		else
		{
			slide_into = spptr->slide_axis().negate();
		}
		
		slide_away = slide_into.negate();
		
		protocols::rigid::RigidBodyTransMoverOP slide_away_from_surface( new protocols::rigid::RigidBodyTransMover( slide_away, pose.num_jump()) );
		slide_away_from_surface->step_size(20);
		
		protocols::docking::FaDockingSlideIntoContactOP slide_into_surface( new protocols::docking::FaDockingSlideIntoContact( pose.num_jump(), slide_into) );

		
		slide_away_from_surface->apply( pose );
		core::Real score_after_slide_away = (*sfxn)(pose);
		TR<<"score after slide away: "<<score_after_slide_away<<std::endl;
		pose.dump_pdb("/Users/mpacella/Rosetta_Surface_Test/after_slide_away.pdb");
		
		//coords before
		//core::Vector starting_coords = pose.residue(1).xyz("CA");
		
		//slide_away_from_surface->apply(pose);
		
		//coords after
		//core::Vector ending_coords = pose.residue(1).xyz("CA");
		
		//TS_ASSERT((ending_coords-starting_coords).normalized() == spptr->slide_axis().negate())
		
		//pose.dump_pdb("/Users/mpacella/Rosetta_Surface_Test/slide_away_from_surf_test.pdb");
		
		
		slide_into_surface->apply(pose);
		core::Real score_after_slide_into = (*sfxn)(pose);
		TR<<"score after slide into: "<<score_after_slide_into<<std::endl;
		pose.dump_pdb("/Users/mpacella/Rosetta_Surface_Test/after_slide_into_surf.pdb");

		protocols::surface_docking::SurfaceOrientMoverOP surf_orient( new protocols::surface_docking::SurfaceOrientMover() );
		surf_orient->set_surface_parameters(spptr);
		surf_orient->apply(pose);
		
		core::Vector coords_after = pose.residue(pose.total_residue()).xyz("Ca2p");
		
		//TR<<"coords_after-coords_before: "<<coords_after-coords_before<<std::endl;
		TR<<"magnitude: "<<(coords_after-coords_before).magnitude()<<std::endl;
															  
		core::Real score_after = (*sfxn)(pose);
		TR<<"score_after_orient: "<<score_after<<std::endl;
		TS_ASSERT(score_after == score_before)
		pose.dump_pdb("/Users/mpacella/Rosetta_Surface_Test/surf_orient_test_after.pdb");
		


	}


};
