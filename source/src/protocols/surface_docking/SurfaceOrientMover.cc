// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SurfaceOrientMover.cc
/// @author Robin A Thottungal (raugust1@jhu.edu)
/// @author Michael Pacella (mpacella88@gmail.com)


// Unit Headers
#include <protocols/surface_docking/SurfaceOrientMover.hh>
// Package Headers
#include <protocols/surface_docking/SurfaceParameters.hh>
// Project headers
#include <core/pose/Pose.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
//Utility Headers
#include <utility/exit.hh>
//Basic Headers
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/xyzVector.hh>
//C++ Headers
#include <math.h>
#include <numeric/xyz.io.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.SurfaceDocking.SurfaceOrientMover" );

namespace protocols {
namespace surface_docking {

using namespace core;
using namespace protocols::moves;
using namespace protocols::surface_docking;

//constructor
SurfaceOrientMover::SurfaceOrientMover():Mover()
	{
    Mover::type( "SurfaceOrientMover");
	}

//destructor
SurfaceOrientMover::~SurfaceOrientMover() {}

void SurfaceOrientMover::apply(pose::Pose & pose)
	{
	//displacement vector from protein centroid to surface centroid
	Vector ProteinCG, SurfaceCG;
	Size const rb_jump=pose.num_jump(); //last jump connects protein to surf
	protocols::geometry::centroids_by_jump (pose, rb_jump, SurfaceCG, ProteinCG);
	TR<<"Initial Protein Centroid:"<<ProteinCG<<std::endl;
	TR<<"Surface Centroid:"<<SurfaceCG<<std::endl;
	Vector total_displacement = ProteinCG - SurfaceCG;
	Vector recenter_vector = calculate_recenter_vector(total_displacement);
	Real recenter_step_size = recenter_vector.magnitude();
	protocols::rigid::RigidBodyTransMover recenter_trans_mover = protocols::rigid::RigidBodyTransMover( recenter_vector, pose.num_jump());
	TR<<"Recenter vector: "<<recenter_vector<<std::endl;
	TR<<"Recenter vector magnitude: "<<recenter_vector.magnitude()<<std::endl;

	recenter_trans_mover.step_size(recenter_step_size);
	recenter_trans_mover.apply(pose);
	TR<<"SurfaceOrient Complete..."<<std::endl;
	}

//TODO: make sure everywhere surface orient apply is called it is followed by slide into contact since that has changed, add check to insure that the score is unchanged by move


void SurfaceOrientMover::set_surface_parameters(protocols::surface_docking::SurfaceParametersOP surface_parameters)
	{
		surface_parameters_ = surface_parameters;
	}

std::string SurfaceOrientMover::get_name() const
	{
		return "SurfaceOrientMover";
	}

Vector SurfaceOrientMover::calculate_recenter_vector( core::Vector const & total_displacement )
	{
		// Change input vector from 'protein displacement from origin' to 'origin displacement from protein'
		core::Vector neg_disp = total_displacement.negated();

		//component of total_displacement along the AB surface vector/////
		Real AB_displacement_component = neg_disp.dot(surface_parameters_->vecAB().normalized());

		//component of total_displacement along the AC surface vector/////
		Real AC_displacement_component = neg_disp.dot(surface_parameters_->vecAC().normalized());

		//how many unit cells are we off of center from in the AB direction?
		SSize num_AB_unit_cells_displaced = SSize(AB_displacement_component/surface_parameters_->vecAB().magnitude());

		//how many unit cells are we off of center from in the AC direction?
		SSize num_AC_unit_cells_displaced = SSize(AC_displacement_component/surface_parameters_->vecAC().magnitude());

	 	//needed AB translation to recenter
		Vector AB_recenter_vector = num_AB_unit_cells_displaced * surface_parameters_->vecAB();
		//needed AC translation to recenter
		Vector AC_recenter_vector = num_AC_unit_cells_displaced * surface_parameters_->vecAC();

		//net translation needed to recenter
		Vector recenter_vector = AB_recenter_vector + AC_recenter_vector;
		return recenter_vector;
	}


}	//surfaceDockingProtocol

}	//protocol
