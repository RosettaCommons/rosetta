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

/// @file SurfaceOrientMover.cc
/// @author Robin A Thottungal (raugust1@jhu.edu)


// Unit Headers
#include <protocols/surface_docking/SurfaceOrientMover.hh>

// Package Headers

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/geometry/RB_geometry.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>

//Utility Headers
#include <utility/exit.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <algorithm>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.SurfaceDocking.SurfaceOrientMover");

////using core::pose::datacache::CacheableDataType::SURFACE_PARAMS;

namespace protocols {
namespace surface_docking {

using namespace core;
using namespace protocols::moves;
using namespace protocols;
//constructor
SurfaceOrientMover::SurfaceOrientMover():Mover(){
    Mover::type( "SurfaceOrientMover");
	}

//destructor
SurfaceOrientMover::~SurfaceOrientMover() {}

void SurfaceOrientMover::apply(pose::Pose & pose){
    // Making sure Surface Vectors are read in properly
    //assert( pose.data().get_ptr( SURFACE_PARAMS )()!=NULL);
    // Calculate the Centroid of the protein and surface ProteinCG, SurfaceCG
    Vector ProteinCG, SurfaceCG;
    Size const rb_jump=pose.num_jump();
    // Last jump is the one that connect the surface with protein
    TR<<"Number of Jumps:"<<pose.num_jump()<<std::endl;
    protocols::geometry::centroids_by_jump (pose, rb_jump,SurfaceCG,ProteinCG);
    TR<<"Initial Protein Centroid:"<<ProteinCG<<std::endl;
    TR<<"Surface Centroid:"<<SurfaceCG<<std::endl;
    // Reading datacache for surface Vectors
    //scoring::SurfaceParameters & surfaceVectors=
    //       *( static_cast< core::scoring::SurfaceParameters * >
	//					( pose.data().get_ptr( SURFACE_PARAMS )() ));
    surface_docking::SurfaceParametersOP surfaceVectors=
    									new surface_docking::SurfaceParameters();
    //Logic Used: If the AB vector is zero, the the surfaceparameters has
    //not been created yet. This will make sure that the surfaceparamters are
    //only calculated once.
    if (surfaceVectors->vecAB.x()==0.0){
        surfaceVectors->GenerateSurfaceParameters(SurfaceCG);
    }

    // Calculating the point of intersection of proteinCG and the plane
    // described by the surface basis vectors
    Vector surfacePointIntersection=
        					surfaceVectors->PlanePointIntersection(ProteinCG);

    TR<<"SurfacePointIntersection:"<<surfacePointIntersection<<std::endl;

    // Calculating the distance between the proteinCG and the surfacePoint
    // intersection (how far is the proteinCG from the surface)
    Real ProteinZdistance=surfacePointIntersection.distance(ProteinCG);

    Vector AC_Intersect=
        surfaceVectors->CalcAxisIntersect(surfacePointIntersection,SurfaceCG,
								surfaceVectors->vecAC,surfaceVectors->vecAB);

    Vector AB_Intersect=
        surfaceVectors->CalcAxisIntersect(surfacePointIntersection,SurfaceCG,
								surfaceVectors->vecAB,surfaceVectors->vecAC);
    TR<<"AC_Intersect:"<<AC_Intersect<<std::endl;
    TR<<"AB_Intersect:"<<AB_Intersect<<std::endl;

    // Calculating the distance between the SurfaceCG and the point that is
    // parallel to the AB & AC vector

    Real ProjectionDistanceAC=AC_Intersect.distance(surfacePointIntersection);

    Real ProjectionDistanceAB=AB_Intersect.distance(surfacePointIntersection);
    // Identifying the translation vector!!!

    Real ABdistance=surfaceVectors->A.distance(surfaceVectors->B);
    Real ACdistance=surfaceVectors->A.distance(surfaceVectors->C);
    TR<<"AC Distance:"<<ACdistance<<std::endl;
    TR<<"AB Distance:"<<ABdistance<<std::endl;
    Vector Trans_AB_Vec,Trans_AC_Vec;

    if (ProjectionDistanceAB <= ABdistance/2.0
        && ProjectionDistanceAC <= ACdistance/2.0) {
        ProteinCG = ProteinCG;
        Trans_AB_Vec=0;
        Trans_AB_Vec=0;
        }

    Trans_AB_Vec=CalcTransVec(ProjectionDistanceAB,
                                        ABdistance,surfaceVectors->vecAB);
    Trans_AC_Vec=CalcTransVec(ProjectionDistanceAC,
                                        ACdistance,surfaceVectors->vecAC);

    TR<<"Trans AC Vec:"<<Trans_AC_Vec<<std::endl;
    TR<<"Trans AB Vec:"<<Trans_AB_Vec<<std::endl;
    TR<<"ProteinZdistance:"<<ProteinZdistance<<std::endl;

    // Add the vector to the point A to get the new co-ordinates
    Vector NewProteinCG;
    NewProteinCG.x()=Trans_AB_Vec.x()+SurfaceCG.x();
    NewProteinCG.y()=Trans_AC_Vec.y()+SurfaceCG.y();
    NewProteinCG.z()=ProteinZdistance+SurfaceCG.z();
    // Identify the Z coordinate on the surface and moving ProteinZdistance
    // along the normal vector axis
    Real a=surfaceVectors->surfacePlane[1];
    Real b=surfaceVectors->surfacePlane[2];
    Real c=surfaceVectors->surfacePlane[3];
    Real d=surfaceVectors->surfacePlane[4];
    Real ZCord4ProteinOnSurface=-((a*NewProteinCG.x())+(b*NewProteinCG.y())
                        +d)/c;
    NewProteinCG.z()=ZCord4ProteinOnSurface;
    Vector UnitNormal=surfaceVectors->unitsurfaceNormalVec;
    NewProteinCG=NewProteinCG+(ProteinZdistance*UnitNormal);

    TR<<"Tramslated ProteinCG:"<<NewProteinCG<<std::endl;
    // For Testing Purpose
    Vector surfaceNewProteinCGIntersection=
                        surfaceVectors->PlanePointIntersection(NewProteinCG);
    TR<<"Intersection point of Tramslated ProteinCG on the plane:"
                                <<surfaceNewProteinCGIntersection<<std::endl;
    // Calculating the distance between the ProteinCentroid and
    // New_ProteinCentroid
    Real trans_magnitude=ProteinCG.distance(NewProteinCG);
    // Moves trans_magnitude distance along the trans_axis
    RigidBodyTransMoverOP TransMover( new RigidBodyTransMover(pose,rb_jump));
    TransMover->trans_axis(NewProteinCG-ProteinCG);
    TransMover->step_size( trans_magnitude );
    TransMover->apply(pose);

    // Setting the transAxis
    TR<<"Translated ProteinCGIntersection:"<<
                  surfaceNewProteinCGIntersection-NewProteinCG<<std::endl;
    surfaceVectors->slideaxis=surfaceNewProteinCGIntersection-NewProteinCG;
    TR<<"SlideAxis:"<<surfaceVectors->slideaxis<<std::endl;
    }

Vector SurfaceOrientMover::CalcTransVec(Real ProjectionDistance,Real
                                            VectorDistance, Vector Vec){
	Real bb,ibb,floor_ceil_ibb;
	Vector TransVec;
	TransVec=0;
	if ( ProjectionDistance > VectorDistance/2.0 ){
		bb = ProjectionDistance / VectorDistance ;
		ibb = floor(bb);
		floor_ceil_ibb = bb - ibb;
		//if (ibb >= 1.0){
		if (floor_ceil_ibb <= 0.5){
			TransVec=floor_ceil_ibb*Vec; // Scaling the vector
		}
		else{
			TransVec=floor_ceil_ibb*Vec; // Scaling the vector
			}
		}
		//}
	return TransVec;
}

std::string SurfaceOrientMover::get_name() const {
	return "SurfaceOrientMover";
	}



}	//surfaceDockingProtocol

}	//protocol
