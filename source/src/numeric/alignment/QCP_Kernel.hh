#ifndef _numeric_alignment_QCP_Kernal_HPP_
#define _numeric_alignment_QCP_Kernal_HPP_

#include <utility/numbers.hh>

#include <cstdlib>
#include <cmath>

namespace numeric
{
namespace alignment
{

template <class Real>
class QCP_Kernel
{

public:
	QCP_Kernel() {}
	virtual ~QCP_Kernel() {}

	///////////////////////////////////////////////////////////////
	/// \remarks
	/// Removes center of mass from given coordinate array.
	///
	/// \param 	coords [In] Array containing coordinates.
	///
	/// \param 	number_of_atoms [In] Number of atoms of both conformations.
	///
	/// \author fordas@uw.edu
	/// \date 10/10/2013
	///////////////////////////////////////////////////////////////
	static void remove_center_of_mass(Real* coordinates, int number_of_atoms)
	{
		Real center[3];
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;

		for(int n = 0; n < number_of_atoms * 3; n += 3)
		{
			center[0] += coordinates[n + 0];
			center[1] += coordinates[n + 1];
			center[2] += coordinates[n + 2];
		}

		center[0] /= number_of_atoms;
		center[1] /= number_of_atoms;
		center[2] /= number_of_atoms;

		for(int n = 0; n < number_of_atoms * 3; n += 3)
		{
			coordinates[n + 0] -= center[0];
			coordinates[n + 1] -= center[1];
			coordinates[n + 2] -= center[2];
		}
	}

	///////////////////////////////////////////////////////////////
	/// \remarks
	/// Implements the 'inner product' operation of Douglas Theobald QCP superposition method (see : http://theobald.brandeis.edu/qcp/
	/// and "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial."  Acta Crystallogr A 61(4):478-480
	/// for more info).
	///
	/// \param 	A [In/Out] 3x3 matrix for the coordinate inner product.
	///
	/// \param 	coords_a [In] Array containing centered coordinates.
	///
	/// \param 	coords_b [In] Array containing centered coordinates.
	///
	/// \param 	number_of_atoms [In] Number of atoms of both conformations.
	///
	/// \return The E0 parameter (upper bound for max Eigenvalue).
	///
	/// \author victor_gil
	/// \date 05/10/2012
	///////////////////////////////////////////////////////////////
	static Real inner_product(
	  Real* A,
	  Real* first_conformation_coords,
	  Real* second_conformation_coords,
	  int number_of_atoms
	)
	{
		Real x1, x2, y1, y2, z1, z2;
		int i;
		Real G1 = 0.0, G2 = 0.0;

		A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

		int total_number_of_coordinates = 3* number_of_atoms;

		for (i = 0; i < total_number_of_coordinates; i+=3)
		{
			x1 = first_conformation_coords[i];
			y1 = first_conformation_coords[i+1];
			z1 = first_conformation_coords[i+2];

			x2 = second_conformation_coords[i];
			y2 = second_conformation_coords[i+1];
			z2 = second_conformation_coords[i+2];

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

		return (G1 + G2) * 0.5;
	}

	///////////////////////////////////////////////////////////////
	/// \remarks
	///	This function ports the second part of Douglas Theobald's QCP superposition method.
	///
	/// \param 	A [In] 3x3 coordinate inner product matrix, as computed via innerProduct.
	///
	/// \param 	E0 [In] Upper bound for the maximum eigenvalue.
	///
	/// \param 	number_of_atoms [In] Number of atoms of conformations used to generate A and E0.
	//
	/// \param 	rot_matrix [Out] 3x3 Output rotation matrix of superposition, if non-null.
	///
	/// \return Rmsd between source conformations.
	///
	/// \author victor_gil
	/// \date 05/10/2012
	///////////////////////////////////////////////////////////////
	static Real calc_rmsd_Theobald_method(
	  Real *A,
	  Real E0,
	  int number_of_atoms,
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
		for (int i = 0; i < 50; ++i)
		{
			oldg = mxEigenV;
			x2 = mxEigenV*mxEigenV;
			b = (x2 + C[2])*mxEigenV;
			a = b + C[1];
			delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
			mxEigenV -= delta;
			if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
				break;
		}

		if (rot_matrix != NULL )
		{
			Real a11, a12, a13, a14, a21, a22, a23, a24,
			     a31, a32, a33, a34, a41, a42, a43, a44;

			Real a3344_4334, a3244_4234, a3243_4233,
			     a3143_4133, a3144_4134, a3142_4132;

			Real q1, q2, q3, q4, normq;

			Real evecprec = 1e-6;

			Real a2, x2, y2, z2;

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

			if (qsqr < evecprec)
			{
				q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
				q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
				q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
				q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
				qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

				if (qsqr < evecprec)
				{
					Real a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
					Real a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
					Real a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

					q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
					q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
					q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
					q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
					qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

					if (qsqr < evecprec)
					{
						q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
						q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
						q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
						q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
						qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;

						if (qsqr < evecprec)
						{
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
		}

		if (utility::is_nan(mxEigenV))
		{
			return 0.0;
		}
		else
		{
			return sqrt(fabs(2.0 * (E0 - mxEigenV)/number_of_atoms));
		}
	}

	///////////////////////////////////////////////////////////////
	/// \remarks
	/// Wrapping function for Douglas Theobald QCP superposition method to calculate the RMSD for two conformations.
	///
	/// \param 	first_conformation_coords [In] Array containing the coordinates of the reference conformation.
	///
	/// \param 	second_conformation_coords [In] Array containing the coordinates of the conformation to be measured.
	///
	/// \param 	number_of_atoms [In] Number of atoms of both conformations.
	///
	/// \param 	rot_matrix [Out] 3x3 Output rotation matrix of superposition, if non-null.
	///
	/// \return The rmsd between both conformations.
	///
	/// \author victor_gil
	/// \date 05/10/2012
	///////////////////////////////////////////////////////////////
	static Real calc_coordinate_rmsd(
	  Real* coords_a,
	  Real* coords_b,
	  int number_of_atoms,
	  Real* rot_matrix)
	{
		remove_center_of_mass(coords_a, number_of_atoms);
		remove_center_of_mass(coords_b, number_of_atoms);

		Real A[9];
		Real E0 = inner_product(A, coords_a, coords_b, number_of_atoms);
		return calc_rmsd_Theobald_method(A, E0, number_of_atoms, rot_matrix);
	}

	///////////////////////////////////////////////////////////////
	/// \remarks
	/// Wrapping function for Douglas Theobald QCP superposition method to calculate the RMSD for two conformations.
	///
	/// \param 	first_conformation_coords [In] Array containing the coordinates of the reference conformation.
	///
	/// \param 	second_conformation_coords [In] Array containing the coordinates of the conformation to be measured.
	///
	/// \param 	number_of_atoms [In] Number of atoms of both conformations.
	///
	/// \param 	rot_matrix [Out] 3x3 Output rotation matrix of superposition, if non-null.
	///
	/// \return The rmsd between both conformations.
	///
	/// \author victor_gil
	/// \date 05/10/2012
	///////////////////////////////////////////////////////////////
	static Real calc_centered_coordinate_rmsd(
	  Real* coords_a,
	  Real* coords_b,
	  int number_of_atoms,
	  Real* rot_matrix)
	{
		Real A[9];
		Real E0 = inner_product(A, coords_a, coords_b, number_of_atoms);
		return calc_rmsd_Theobald_method(A, E0, number_of_atoms, rot_matrix);
	}

};
}
}

#endif
