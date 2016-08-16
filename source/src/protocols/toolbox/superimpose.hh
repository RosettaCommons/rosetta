// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/Job.hh
/// @brief  header file for ThreadingJob classes, part of August 2008 job distributor as planned at RosettaCon08.  This file is responsible for three ideas: "inner" jobs, "outer" jobs (with which the job distributor works) and job container (currently just typdefed in the .fwd.hh)
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_toolbox_superimpose_hh
#define INCLUDED_protocols_toolbox_superimpose_hh

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>


namespace protocols {
namespace toolbox {

typedef  numeric::xyzMatrix< core::Real > Matrix;
typedef  numeric::xyzVector< core::Real > Vector;

void CA_superimpose( ObjexxFCL::FArray1_double const& weights, core::pose::Pose const& ref_pose, core::pose::Pose& fit_pose );
void CA_superimpose( core::pose::Pose const& ref_pose, core::pose::Pose& fit_pose );
/*
void superimpose(
core::Size natoms,
ObjexxFCL::FArray1_double& weights,
ObjexxFCL::FArray2_double& ref_coords,
ObjexxFCL::FArray2_double& coords,
Matrix &R //returns rotation matrix
);

void superimpose(
core::Size natoms,
ObjexxFCL::FArray1_double& weights,
ObjexxFCL::FArray2_double& ref_coords,
ObjexxFCL::FArray2_double& coords
); */

/* @brief Calculates transform need to superimpose init_coords onto ref_coords.
*
* Returned transform types:
*  to_init_center - Transform placing init_coords center at origin. Ie. -(init_coords center point).
*  to_fit_center -  Transform placing ref_coords center at origin. Ie. -(ref_coords center point).
*  rotation - Rotation matrix applied about center move init to fit.
*
*  Ex. To apply superposition transform:
*
*   fit_coordiate_value = to_fit_center + rotation * (init_coordinate_value - to_init_center)
*/
void
superposition_transform(
	utility::vector1< numeric::xyzVector< core::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< core::Real > > & ref_coords,
	Matrix & rotation,
	Vector & to_init_center,
	Vector & to_fit_center);

/* @brief Calculates transform need to superimpose init_coords onto ref_coords using the given coordinate weights.
*
* Returned transform types:
*  to_init_center - Transform placing init_coords center at origin. Ie. -(init_coords center point).
*  to_fit_center -  Transform placing ref_coords center at origin. Ie. -(ref_coords center point).
*  rotation - Rotation matrix applied about center move init to fit.
*
*  Ex. To apply superposition transform:
*
*   fit_coordiate_value = to_fit_center + rotation * (init_coordinate_value - to_init_center)
*/
void
superposition_transform(
	utility::vector1< numeric::xyzVector< core::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< core::Real > > & ref_coords,
	utility::vector1< core::Real >   & coord_weights,
	Matrix & rotation,
	Vector & to_init_center,
	Vector & to_fit_center);

/* @brief Applies given superposition transform series to pose jump.
*
* Applies superposition to jump, updating all downstream pose components.
*/
void apply_superposition_transform_to_jump(
	core::pose::Pose & pose,
	core::Size jump_id,
	Matrix rotation,
	Vector to_init_center,
	Vector to_fit_center);

/* @brief Applies given superposition transform series to pose.
*
* Applies superposition transform to all atoms within pose.
*/
void apply_superposition_transform(
	core::pose::Pose & pose,
	Matrix rotation,
	Vector to_init_center,
	Vector to_fit_center);

/* @brief Converts a vector1-of-xyzVectors into FArray2D. */
template <typename T>
void vector_vector_to_FArray2(
	utility::vector1< numeric::xyzVector< T > > & from,
	ObjexxFCL::FArray2D< T > & to);

/* @brief Calculates superposition transform from coords to ref_coords.
*
* Modifies ref_coords and coords, moving coords into superposition.
*/
void superposition_transform(
	core::Size natoms,
	ObjexxFCL::FArray1_double const& weights,
	ObjexxFCL::FArray2_double& ref_coords,
	ObjexxFCL::FArray2_double& coords,
	Matrix &R,
	Vector &toCenter,
	Vector &toFitCenter);

void fit_centered_coords(
	core::Size natoms,
	ObjexxFCL::FArray1_double const& weights,
	ObjexxFCL::FArray2_double const& ref_coords,
	ObjexxFCL::FArray2_double& coords,
	Matrix &R
);

void reset_x(
	core::Size n,
	ObjexxFCL::FArray2_double& x,
	ObjexxFCL::FArray1_double const& wts,
	ObjexxFCL::FArray1_double& transvec
);

/// @brief write a CA ALA pdb
void dump_as_pdb(
	std::string filename,
	core::Size n,
	ObjexxFCL::FArray2_double& coords,
	ObjexxFCL::FArray1D_double transvec = ObjexxFCL::FArray1D_double( 3, 0.0 )
);

void fill_CA_coords( core::pose::Pose const& pose, core::Size natoms, ObjexxFCL::FArray2_double& coords );
void fill_CA_coords( core::pose::Pose const& pose, ObjexxFCL::FArray2_double& coords );

}
}

#endif
