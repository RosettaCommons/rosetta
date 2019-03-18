// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// (c) Components copyright (c) 2009-2016 Pu Liu and Douglas L. Theobald
// (c) All rights reserved.
//
// (c) Redistribution and use in source and binary forms, with or without modification, are permitted
// (c) provided that the following conditions are met:
//
// (c) * Redistributions of source code must retain the above copyright notice, this list of
// (c)   conditions and the following disclaimer.
// (c) * Redistributions in binary form must reproduce the above copyright notice, this list
// (c)   of conditions and the following disclaimer in the documentation and/or other materials
// (c)   provided with the distribution.
//
// (c) THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// (c) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// (c) LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// (c) PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// (c) HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// (c) SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// (c) LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// (c) DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// (c) THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (c) (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// (c) OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _numeric_alignment_QCPKernel_HPP_
#define _numeric_alignment_QCPKernel_HPP_

#include <cstdlib>
#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace numeric
{
namespace alignment
{

template <class Real>
class QCPKernel
{  public:
	typedef Eigen::Matrix< Real, 3, 1 >  Point;
	typedef Eigen::Matrix< Real, 3, Eigen::Dynamic >  Coords;
	typedef Eigen::Map< Coords >  CoordMap;

	/*
	* Implements the 'inner product' operation of Douglas Theobald QCP superposition method (see : http://theobald.brandeis.edu/qcp/
	* and "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial."  Acta Crystallogr A 61(4):478-480
	* for more info).
	*
	* @param  A 3x3 matrix buffer coordinate inner product.
	* @param  E0 Upper bound for max Eigenvalue.
	*
	* @param  first_coords Coordinate matrix for alignment target coordinates.
	* @param  first_coords_center_of_mass Center of mass for first_coords.
	*
	* @param  second_coords Coordinate matrix for alignment source coordinates.
	* @param  second_coords_center_of_mass Center of mass for second_coords.
	*/
	template <typename DerivedA, typename DerivedB>
	static void inner_product(
		Real* A,
		Real& E0,
		const Eigen::DenseBase<DerivedA>& first_coords,
		const Point & first_coords_center_of_mass,
		const Eigen::DenseBase<DerivedB>& second_coords,
		const Point & second_coords_center_of_mass
	)
	{
		Real x1, x2, y1, y2, z1, z2;
		size_t i;
		Real G1 = 0.0, G2 = 0.0;

		A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

		size_t total_number_of_coordinates = first_coords.cols();

		for ( i = 0; i < total_number_of_coordinates; i++ ) {
			x1 = first_coords.col(i)[0] - first_coords_center_of_mass[0];
			y1 = first_coords.col(i)[1] - first_coords_center_of_mass[1];
			z1 = first_coords.col(i)[2] - first_coords_center_of_mass[2];

			x2 = second_coords.col(i)[0] - second_coords_center_of_mass[0];
			y2 = second_coords.col(i)[1] - second_coords_center_of_mass[1];
			z2 = second_coords.col(i)[2] - second_coords_center_of_mass[2];

			G1 += x1 * x1 + y1 * y1 + z1 * z1;
			G2 += x2 * x2 + y2 * y2 + z2 * z2;

			A[0] +=  (x1 * x2);
			A[1] +=  (x1 * y2);
			A[2] +=  (x1 * z2);

			A[3] +=  (y1 * x2);
			A[4] +=  (y1 * y2);
			A[5] +=  (y1 * z2);

			A[6] +=  (z1 * x2);
			A[7] +=  (z1 * y2);
			A[8] +=  (z1 * z2);
		}

		E0 = (G1 + G2) * .5;
	}

	/*
	* Second component of Douglas Theobald's QCP superposition method.
	*
	* @param  A 3x3 coordinate inner product matrix, as computed via inner_product.
	* @param  E0 Upper bound for the maximum eigenvalue, as computed via inner_product.
	* @param  number_of_atoms Number of atoms of conformations used to generate A and E0.
	*
	* @param  rot_matrix 3x3 Output rotation matrix of superposition, if non-null.
	*
	* @return RMSD between conformations.
	*/
	static Real calc_rmsd_Theobald_method(
		Real *A,
		Real E0,
		size_t number_of_atoms,
		Real* rot_matrix
	)
	{
		Real Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
		Real Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
			SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
			SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
			SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
		Real C[4];
		Real mxEigenV;
		Real b, a, delta, x2;
		Real oldg = 0.0;
		Real evalprec = 1e-11;

		Sxx = A[0];
		Sxy = A[1];
		Sxz = A[2];
		Syx = A[3];
		Syy = A[4];
		Syz = A[5];
		Szx = A[6];
		Szy = A[7];
		Szz = A[8];

		Sxx2 = Sxx * Sxx;
		Syy2 = Syy * Syy;
		Szz2 = Szz * Szz;

		Sxy2 = Sxy * Sxy;
		Syz2 = Syz * Syz;
		Sxz2 = Sxz * Sxz;

		Syx2 = Syx * Syx;
		Szy2 = Szy * Szy;
		Szx2 = Szx * Szx;

		SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
		Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

		C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
		C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

		SxzpSzx = Sxz + Szx;
		SyzpSzy = Syz + Szy;
		SxypSyx = Sxy + Syx;
		SyzmSzy = Syz - Szy;
		SxzmSzx = Sxz - Szx;
		SxymSyx = Sxy - Syx;
		SxxpSyy = Sxx + Syy;
		SxxmSyy = Sxx - Syy;
		Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

		C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
			+ (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
			+ (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
			+ (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
			+ (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
			+ (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

		mxEigenV = E0;
		for ( size_t i = 0; i < 50; ++i ) {
			oldg = mxEigenV;
			x2 = mxEigenV*mxEigenV;
			b = (x2 + C[2])*mxEigenV;
			a = b + C[1];
			delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
			mxEigenV -= delta;
			if ( fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV) ) {
				break;
			}
		}

		if ( rot_matrix != NULL ) {
			Real a11, a12, a13, a14, a21, a22, a23, a24,
				a31, a32, a33, a34, a41, a42, a43, a44;

			Real a3344_4334, a3244_4234, a3243_4233,
				a3143_4133, a3144_4134, a3142_4132;

			Real q1, q2, q3, q4, normq;

			Real evecprec = 1e-6;

			Real a2, y2, z2;

			Real xy, az, zx, ay, yz, ax;

			a11 = SxxpSyy + Szz - mxEigenV;
			a12 = SyzmSzy;
			a13 = -SxzmSzx;
			a14 = SxymSyx;
			a21 = SyzmSzy;
			a22 = SxxmSyy - Szz-mxEigenV;
			a23 = SxypSyx;
			a24= SxzpSzx;
			a31 = a13;
			a32 = a23;
			a33 = Syy-Sxx-Szz - mxEigenV;
			a34 = SyzpSzy;
			a41 = a14;
			a42 = a24;
			a43 = a34;
			a44 = Szz - SxxpSyy - mxEigenV;
			a3344_4334 = a33 * a44 - a43 * a34;
			a3244_4234 = a32 * a44-a42*a34;
			a3243_4233 = a32 * a43 - a42 * a33;
			a3143_4133 = a31 * a43-a41*a33;
			a3144_4134 = a31 * a44 - a41 * a34;
			a3142_4132 = a31 * a42-a41*a32;
			q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
			q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
			q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
			q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

			Real qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

			if ( qsqr < evecprec ) {
				q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
				q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
				q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
				q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
				qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

				if ( qsqr < evecprec ) {
					Real a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
					Real a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
					Real a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

					q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
					q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
					q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
					q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
					qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

					if ( qsqr < evecprec ) {
						q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
						q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
						q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
						q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
						qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;

						if ( qsqr < evecprec ) {
							// if qsqr is still too small, return the identity matrix.
							rot_matrix[0] = rot_matrix[4] = rot_matrix[8] = 1.0;
							rot_matrix[1] = rot_matrix[2] = rot_matrix[3] = rot_matrix[5] = rot_matrix[6] = rot_matrix[7] = 0.0;
						}
					}
				}
			}

			normq = sqrt(qsqr);
			q1 /= normq;
			q2 /= normq;
			q3 /= normq;
			q4 /= normq;

			a2 = q1 * q1;
			x2 = q2 * q2;
			y2 = q3 * q3;
			z2 = q4 * q4;

			xy = q2 * q3;
			az = q1 * q4;
			zx = q4 * q2;
			ay = q1 * q3;
			yz = q3 * q4;
			ax = q1 * q2;

			rot_matrix[0] = a2 + x2 - y2 - z2;
			rot_matrix[1] = 2 * (xy + az);
			rot_matrix[2] = 2 * (zx - ay);
			rot_matrix[3] = 2 * (xy - az);
			rot_matrix[4] = a2 - x2 + y2 - z2;
			rot_matrix[5] = 2 * (yz + ax);
			rot_matrix[6] = 2 * (zx + ay);
			rot_matrix[7] = 2 * (yz - ax);
			rot_matrix[8] = a2 - x2 - y2 + z2;
		}// if rot_matrix != NULL

		if ( std::isnan(mxEigenV) ) {
			return 0.0;
		} else {
			return sqrt(fabs(2.0 * (E0 - mxEigenV)/number_of_atoms));
		}
	}

	/*
	* Calculate rmsd between given coordinate sets via QCP superposition method.
	*
	* @param  first_coords                Shape (3,n) matrix type containing target coordinates.
	* @param  second_coords                Shape (3,n) matrix type containing target coordinates.
	*
	* @return RMSD between tconformations.
	*/
	template <typename DerivedA, typename DerivedB>
	static Real calc_coordinate_rmsd(
		const Eigen::DenseBase<DerivedA>& first_coords,
		const Eigen::DenseBase<DerivedB>& second_coords
	)
	{
		Point first_coords_center_of_mass = first_coords.rowwise().sum() / first_coords.cols();
		Point second_coords_center_of_mass = second_coords.rowwise().sum() / second_coords.cols();

		return calc_coordinate_rmsd(
			first_coords, first_coords_center_of_mass, second_coords, second_coords_center_of_mass);
	}

	/*
	* Calculate rmsd between given coordinate sets via QCP superposition method.
	*
	* @param  first_coords                Shape (3,n) matrix type containing target coordinates.
	* @param  first_coords_center_of_mass Shape (3,1) center of mass for first_coords.
	*
	* @param  second_coords                Shape (3,n) matrix type containing target coordinates.
	* @param  second_coords_center_of_mass Shape (3,1) center of mass for second_coords.
	*
	* @return RMSD between tconformations.
	*/
	template <typename DerivedA, typename DerivedB>
	static Real calc_coordinate_rmsd(
		const Eigen::DenseBase<DerivedA>& first_coords,
		const Point & first_coords_center_of_mass,
		const Eigen::DenseBase<DerivedB>& second_coords,
		const Point & second_coords_center_of_mass
	)
	{
		Real A[9];
		Real E0;

		inner_product(A, E0, first_coords, first_coords_center_of_mass, second_coords, second_coords_center_of_mass);
		return calc_rmsd_Theobald_method(A, E0, first_coords.cols(), NULL);
	}

	/*
	* Calculate rmsd and superposition transform superimposing src_coords onto onto_coords.
	*
	* @param  src_coords                Shape (3,n) matrix type containing src coordinates.
	* @param  src_coords_center_of_mass Shape (3,1) center of mass for src_coords.
	*
	* @param  onto_coords                Shape (3,n) matrix type containing onto coordinates.
	* @param  onto_coords_center_of_mass Shape (3,1) center of mass for onto_coords.
	*
	* @param superposition_transform      Output superposition transform.
	*
	* @return RMSD between conformations.
	*/
	template <typename DerivedA, typename DerivedB>
	static Real calc_coordinate_superposition(
		const Eigen::DenseBase<DerivedA>& src_coords,
		const Point & src_center_of_mass,
		const Eigen::DenseBase<DerivedB>& onto_coords,
		const Point & onto_center_of_mass,
		Eigen::Transform<Real, 3, Eigen::Affine> & superposition_transform
	)
	{
		Real A[9];
		Real E0;
		Real superposition_rotation[9];
		Eigen::Map< Eigen::Matrix<Real, 3, 3, Eigen::RowMajor> >  rotation_matrix(superposition_rotation);
		Eigen::Translation<Real, 3> src_center(src_center_of_mass);
		Eigen::Translation<Real, 3> onto_center(onto_center_of_mass);

		inner_product(
			A, E0,
			onto_coords, onto_center_of_mass,
			src_coords, src_center_of_mass
		);
		Real rmsd = calc_rmsd_Theobald_method(A, E0, onto_coords.cols(), superposition_rotation);

		superposition_transform = onto_center * rotation_matrix * src_center.inverse();
		return rmsd;
	}

	/*
	* Calculate rmsd and superposition transform superimposing src_coords onto onto_coords.
	*
	* @param  src_coords                Shape (3,n) matrix type containing src coordinates.
	*
	* @param  onto_coords                Shape (3,n) matrix type containing onto coordinates.
	*
	* @param superposition_transform      Output superposition transform.
	*
	* @return RMSD between conformations.
	*/
	template <typename DerivedA, typename DerivedB>
	static Real calc_coordinate_superposition(
		const Eigen::DenseBase<DerivedA>& src_coords,
		const Eigen::DenseBase<DerivedB>& onto_coords,
		Eigen::Transform<Real, 3, Eigen::Affine> & superposition_transform
	)
	{
		Point src_center_of_mass = src_coords.rowwise().sum() / src_coords.cols();
		Point onto_center_of_mass = onto_coords.rowwise().sum() / onto_coords.cols();

		return calc_coordinate_superposition(
			src_coords, src_center_of_mass,
			onto_coords, onto_center_of_mass,
			superposition_transform
		);
	}


};
}
}

#endif
