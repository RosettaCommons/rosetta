// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/LoopGrid.hh
/// @brief  generate accessible grids for loop
/// @author Yuan Liu (wendao@u.washington.edu)

#ifndef INCLUDED_devel_loop_grids_LoopGrid_hh
#define INCLUDED_devel_loop_grids_LoopGrid_hh


// Unit headers
#include <core/pose/Pose.hh>
#include <protocols/loops/Loop.hh>

//#include <protocols/match/LoopGrid.fwd.hh>
#include <devel/loop_grids/LoopGrid.fwd.hh>

#include <protocols/match/BumpGrid.fwd.hh>
#include <utility/vector1.hh>

#include <numeric/geometry/BoundingBox.fwd.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>


namespace protocols {
namespace match {

class LoopGrid : public utility::pointer::ReferenceCount
{
public:
	typedef numeric::geometry::BoundingBox< core::Vector > BoundingBox;
	typedef numeric::geometry::hashing::Bin3D Bin3D;
	typedef core::Vector Vector;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef protocols::loops::Loop Loop;
	typedef protocols::loops::LoopOP LoopOP;
	typedef protocols::loops::LoopCOP LoopCOP;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::pose::PoseCOP PoseCOP;
	typedef core::conformation::Residue Residue;
	typedef utility::fixedsizearray1< Vector, 8 > CornerPoints;
	typedef utility::fixedsizearray1< Real, 2 > xyVector;

public:
	LoopGrid( Pose const &pose, Loop const &loop );

	//if the point or all CAs are in the grid, return true
	//else return false, which means clash
	bool occupied( Vector const & xyz ) const;
	bool occupied( Bin3D const & bin ) const;
	//the whole protein
	bool occupied( Pose const &pose ) const;
	//part of the protein, contain loop
	bool occupied( Pose const &pose, Size start ) const;
	//part of the loop
	bool occupied( Pose const &pose, Size start, Size stop ) const;
	//residue
	bool occupied( Residue const &res, Size nres ) const;

	virtual ~LoopGrid();

private:
	utility::vector1< Bool3DGridOP > grids_; //grids for each res of the loop 2 ~ n-1
	Real bin_width_; //bin width
	Pose pose_; //bump scaffold, exclude the loop
	Loop loop_; //loop
	static const Size MAX_LENTH = 60; //normally, we can not deal with len>12, but i statistic to 60
	static const Real loop_length_cutoff[MAX_LENTH+1]; //length cutoff that contain 99% loop

	/// @brief create the bump grid for specify res, the two centers are the roots of the loop
	/// @brief Radius is value in loop_length_cutoff[], the ith res: loop_left=i loop_right=n-i+1
	void create_grid_and_by_two_sphere( Bool3DGridOP &gridop,
		Vector center1, Real r1, Vector center2, Real r2 );

	/// @brief output the grid of ndx th res, just for debug
	void write_grids_file(const char *fn, Size ndx);

	//CornerPoints overlap_boundary_2d_spheres( xyVector center1, Real r1, xyVector center2, Real r2 );
};//class

}//namespace match
}//namespace protocols

#endif

