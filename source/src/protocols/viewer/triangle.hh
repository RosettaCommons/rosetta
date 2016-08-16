// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_viewer_triangle_hh
#define INCLUDED_protocols_viewer_triangle_hh

#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace viewer {


// 88 bytes
class triangle {
public:
	triangle() {
		depth_ = 0;
	}

	triangle( const numeric::xyzVector<float> vertices[], const numeric::xyzVector<float> normals[] ) {
		for ( int i=0; i<3; ++i ) {
			for ( int j=0; j<3; ++j ) {
				vertices_[i][j] = (float)vertices[i][j];
				normals_[i][j] = (float)normals[i][j];
			}
		}
		center_ = 0.33333f * (vertices_[0]+vertices_[1]+vertices_[2]);
		depth_ = 0;

	}

	triangle( const numeric::xyzVector<double> vertices[], const numeric::xyzVector<double> normals[] ) {
		for ( int i=0; i<3; ++i ) {
			for ( int j=0; j<3; ++j ) {
				vertices_[i][j] = (float)vertices[i][j];
				normals_[i][j] = (float)normals[i][j];
			}
		}
		center_ = 0.33333f * (vertices_[0]+vertices_[1]+vertices_[2]);
		depth_ = 0;
	}

	// ASSUMES zdir IS NORMALIZED!!!!!!!
	void update( const numeric::xyzVector_float &zdir ) {
		//depth_ = numeric::dot( center_, zdir );
		depth_ = center_.x() * zdir.x() + center_.y()*zdir.y() + center_.z()*zdir.z();
	}

	float depth_;
	numeric::xyzVector_float center_;
	numeric::xyzVector_float vertices_[3], normals_[3];
};

}
}

#endif
