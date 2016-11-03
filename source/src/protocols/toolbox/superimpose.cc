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

#if (defined _WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

#include <protocols/toolbox/superimpose.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

#include <core/types.hh>
#include <ObjexxFCL/FArray2D.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <numeric/xyz.functions.hh>
//// C++ headers
#include <string>
#include <iostream>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <utility/vector1.hh>
#include <numeric/model_quality/RmsData.hh>


namespace protocols {
namespace toolbox {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer tr( "protocols.evaluation.PCA", basic::t_info );

using namespace core;
using namespace numeric::model_quality; //for rms functions

typedef core::Real matrix[3][3];
typedef core::Real rvec[3];

template <typename T>
void
vector_vector_to_FArray2(
	utility::vector1< numeric::xyzVector< T > > & from,
	ObjexxFCL::FArray2D< T > & to)
{
	core::Size count = from.size();

	/*
	runtime_assert(to.I1() == 3);
	runtime_assert(to.I2() == count);
	*/

	for ( core::Size i = 1; i <= count; i++ ) {
		for ( core::Size j = 1; j <= 3; j++ ) {
			to(j,i) = from[i](j);
		}
	}
}

void
superposition_transform(
	utility::vector1< numeric::xyzVector< core::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< core::Real > > & ref_coords,
	numeric::xyzMatrix< core::Real > & rotation,
	numeric::xyzVector< core::Real > & to_init_center,
	numeric::xyzVector< core::Real > & to_fit_center)
{
	utility::vector1< core::Real > dummy_weights(0);

	superposition_transform(
		init_coords,
		ref_coords,
		dummy_weights,
		rotation,
		to_init_center,
		to_fit_center);
}

void
superposition_transform(
	utility::vector1< numeric::xyzVector< core::Real > > & init_coords,
	utility::vector1< numeric::xyzVector< core::Real > > & ref_coords,
	utility::vector1< core::Real > & coord_weights,
	numeric::xyzMatrix< core::Real > & rotation,
	numeric::xyzVector< core::Real > & to_init_center,
	numeric::xyzVector< core::Real > & to_fit_center)
{
	core::Size count = init_coords.size();
	runtime_assert(count == ref_coords.size());

	ObjexxFCL::FArray2D< numeric::Real > init_fa( 3, init_coords.size());
	ObjexxFCL::FArray2D< numeric::Real > ref_fa( 3, ref_coords.size());

	vector_vector_to_FArray2(init_coords, init_fa);
	vector_vector_to_FArray2(ref_coords, ref_fa);

	ObjexxFCL::FArray1D< numeric::Real > weights_fa(count, 1);
	if ( coord_weights.size() != 0 ) {
		runtime_assert(count == coord_weights.size());
		for ( core::Size i = 1; i <= count; i++ ) {
			weights_fa[i] =  coord_weights[i];
		}
	}

	protocols::toolbox::superposition_transform(
		count,
		weights_fa,
		ref_fa,
		init_fa,
		rotation,
		to_init_center,
		to_fit_center);
}

void apply_superposition_transform_to_jump(
	core::pose::Pose & pose,
	core::Size jump_id,
	Matrix rotation,
	Vector to_init_center,
	Vector to_fit_center)
{
	using namespace core::kinematics;

	Stub upstream_stub = pose.conformation().upstream_jump_stub(jump_id);
	Jump target_jump = pose.jump(jump_id);

	target_jump.translation_along_axis(
		upstream_stub,
		to_init_center.normalized(),
		to_init_center.length());

	target_jump.rotation_by_matrix(
		upstream_stub,
		Vector(0, 0, 0),
		rotation);

	target_jump.translation_along_axis(
		upstream_stub,
		-to_fit_center.normalized(),
		to_fit_center.length());

	pose.set_jump(jump_id, target_jump);
}

void apply_superposition_transform(
	core::pose::Pose & pose,
	Matrix rotation,
	Vector to_init_center,
	Vector to_fit_center)
{
	for ( Size r = 1; r <= pose.size(); r++ ) {
		for ( Size a = 1; a <= pose.residue_type(r).natoms(); a++ ) {
			core::id::AtomID id(a, r);
			pose.set_xyz( id, ( rotation * ( pose.xyz(id) + to_init_center) ) - to_fit_center );
		}
	}
}

void fill_CA_coords( core::pose::Pose const& pose, FArray2_double& coords ) {  // fill coords
	fill_CA_coords( pose, pose.size(), coords );
}

//coords needs to be a 3 x total_residues array
void fill_CA_coords( core::pose::Pose const& pose, Size n_atoms, FArray2_double& coords ) {  // fill coords
	for ( core::Size i = 1; i <= n_atoms; i++ ) {
		id::NamedAtomID idCA( "CA", i );
		PointPosition const& xyz = pose.xyz( idCA );
		for ( core::Size d = 1; d<=3; ++d ) {
			coords( d, i ) = xyz[ d-1 ];
		}
	}
}

void CA_superimpose( FArray1_double const& weights, core::pose::Pose const&  ref_pose, core::pose::Pose& fit_pose ) {
	//count residues with CA
	Size nres=0;
	for ( core::Size i = 1; i <= fit_pose.size(); i++ ) {
		if ( !fit_pose.residue_type( i ).is_protein() ) break;
		++nres;
	}
	Size const natoms( nres );

	FArray2D_double coords( 3, natoms, 0.0 );
	FArray2D_double ref_coords( 3, natoms, 0.0 );
	fill_CA_coords( ref_pose, natoms, ref_coords );
	fill_CA_coords( fit_pose, natoms, coords );

	FArray1D_double transvec( 3 );
	FArray1D_double ref_transvec( 3 );
	reset_x( natoms, coords, weights, transvec );
	reset_x( natoms, ref_coords, weights, ref_transvec );

	Matrix R;
	fit_centered_coords( natoms, weights, ref_coords, coords, R );

	//now apply transvecs and R to pose, via root!
	//   //ASSUME FOR NOW root is first N
	//   id::NamedAtomID idROOT( "N", 1 );
	//   Vector root_xyz( fit_pose.xyz( idROOT ) );
	//   root_xyz += ;// - Vector( ref_transvec( 1 ), ref_transvec( 2 ), ref_transvec( 3 ));
	//   fit_pose.set_xyz( idROOT, root_xyz );

	Vector toCenter( transvec( 1 ), transvec( 2 ), transvec( 3 ));
	Vector toFitCenter( ref_transvec( 1 ), ref_transvec( 2 ), ref_transvec( 3 ));
	{ // translate xx2 by COM and fill in the new ref_pose coordinates
		Size atomno(0);
		Vector x2;
		for ( Size i=1; i<= fit_pose.size(); ++i ) {
			for ( Size j=1; j<= fit_pose.residue_type(i).natoms(); ++j ) { // use residue_type to prevent internal coord update
				++atomno;
				fit_pose.set_xyz( id::AtomID( j,i), R * ( fit_pose.xyz( id::AtomID( j, i) ) + toCenter ) - toFitCenter );
			}
		}
	}

	if ( tr.Trace.visible() ) {
		FArray2D_double pose_coords( 3, natoms );
		fill_CA_coords( fit_pose, natoms, pose_coords );
		for ( Size pos = 1; pos <= natoms; ++pos ) {
			tr.Trace << " fit_coords vs pose_coords: ";
			for ( Size d = 1; d <= 3; ++d )  tr.Trace << coords( d, pos ) << " ";
			for ( Size d = 1; d <= 3; ++d )  tr.Trace << pose_coords( d, pos ) << " ";
			tr.Trace << std::endl;
		}
	}
}

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
	Vector &toFitCenter)
{
	// Move ref and coords to place center of mass at origin.
	// Save resulting transforms.
	FArray1D_double ref_transvec( 3 );
	reset_x( natoms, ref_coords, weights, ref_transvec );
	toFitCenter = Vector( ref_transvec( 1 ), ref_transvec( 2 ), ref_transvec( 3 ));

	FArray1D_double transvec( 3 );
	reset_x( natoms, coords, weights, transvec );
	toCenter = Vector( transvec( 1 ), transvec( 2 ), transvec( 3 ));

	// Fit centered coords, updates cords array with fit location
	fit_centered_coords( natoms, weights, ref_coords, coords, R );

	// Move superimposed coords into reference location
	for ( core::Size i = 1; i <= natoms; i++ ) {
		for ( core::Size d = 1; d<=3; ++d ) {
			ref_coords( d, i ) = ref_coords( d, i ) - ref_transvec(d);
			coords( d, i ) = coords( d, i ) - ref_transvec(d);
		}
	}
}


/// @brief compute projections for given pose
void superimpose( Size natoms, FArray1_double const& weights, FArray2_double& ref_coords, FArray2_double& coords ) {
	//fill Farray for fit

	FArray1D_double transvec( 3 );
	reset_x( natoms, ref_coords, weights, transvec );
	reset_x( natoms, coords, weights, transvec );
	Matrix R; //to return rotation matrix... but we don't care.
	fit_centered_coords( natoms, weights, ref_coords, coords, R );
}

/// @brief compute projections for given pose
void superimpose( Size natoms, FArray1_double const& weights, FArray2_double& ref_coords, FArray2_double& coords, Matrix &R ) {
	//fill Farray for fit
	FArray1D_double transvec( 3 );
	reset_x( natoms, ref_coords, weights, transvec );
	reset_x( natoms, coords, weights, transvec );
	fit_centered_coords( natoms, weights, ref_coords, coords, R );
}


void CA_superimpose( core::pose::Pose const&  ref_pose, core::pose::Pose& fit_pose ) {
	FArray1D_double weights( ref_pose.size(), 1.0 );
	CA_superimpose( weights, ref_pose, fit_pose );
}


void calc_fit_R(int natoms, Real const* weights, rvec const* xp,rvec const*x, matrix R );
void jacobi(double a[6][6],double d[],double v[6][6],int *nrot);

/// @brief A function (not a macro) that will not print a square matrix to tr.Debug
template< class T > void dump_matrix( Size, T const &, basic::Tracer & ) {}

/// @brief A function (not a macro) that will print a square matrix to tr.Debug
template< class T > void dump_matrix_no( Size nr, T const & a, basic::Tracer & tr)
{
	Size i,k;
	for ( i =0 ; i<nr; ++i ) {
		for ( k =0 ; k<nr; ++k ) {
			tr.Debug << a[i][k] << " ";
		}
		tr.Debug << "\n";
	}
}

/// some low-level helper routines

#define DIM 3
void rotate_vec(int natoms,rvec *x,matrix R)
{
	int j,r,c,m;
	rvec x_old;

	/*rotate X*/
	for ( j=0; j<natoms; j++ ) {
		for ( m=0; m<DIM; m++ ) {
			x_old[m]=x[j][m];
		}
		for ( r=0; r<DIM; r++ ) {
			x[j][r]=0;
			for ( c=0; c<DIM; c++ ) {
				x[j][r]+=R[r][c]*x_old[c];
			}
		}
	}
}

void add_vec( int natoms, rvec *x, rvec transvec ) {
	for ( int i=0; i<natoms; i++ ) {
		for ( int j = 0; j< DIM; j++ ) {
			x[i][j]+=transvec[j];
		}
	}
}

void reset_x( Size n, FArray2_double& x, FArray1_double const& wts, FArray1_double& transvec ) {
	Size const dim( 3 );
	Real mass( 0.0 );
	transvec = FArray1D_double( 3, 0.0 );

	for ( Size j = 1; j <= n; ++j ) {
		mass += wts( j );
		// align center of mass to origin
		for ( Size d = 1; d <= dim; ++d ) {
			transvec( d ) += x( d, j )*wts( j );
		}
	}
	for ( Size d = 1; d <= dim; ++d ) {
		transvec( d ) = -transvec( d )/mass;
	}
	for ( Size j = 1; j <= n; ++j ) {
		for ( Size d = 1; d<= dim; ++d ) {
			x( d, j ) += transvec( d );
		}
	}
}


void dump_as_pdb( std::string filename, Size n, FArray2_double& x,  FArray1D_double transvec ) {
	utility::io::ozstream out( filename );
	for ( Size i = 1; i <= n; ++i ) {
		char buf[400];
		sprintf(buf, "ATOM  %5d  %-4s%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f", (int) i, "CA", "ALA", (int) i,
			x( 1, i ) + transvec( 1 ),
			x( 2, i ) + transvec( 2 ),
			x( 3, i ) + transvec( 3 ), 1.0, 1.0 );
		out << buf << std::endl;
	}

}


void fit_centered_coords( Size natoms, FArray1_double const& weights, FArray2_double const& ref_coords, FArray2_double& coords, Matrix &R ) {
	rvec* xgmx;
	rvec* xrefgmx;
	Real* weights_gmx;

	xgmx = new rvec[ natoms ];
	xrefgmx = new rvec[ natoms ];
	weights_gmx = new Real[ natoms ];
	matrix Rot;

	//transfer into C-style arrays
	for ( Size i = 1; i<=natoms ; i++ ) {
		weights_gmx[i-1] = weights( i );
		for ( Size d = 1; d<=3; d++ ) {
			xgmx[i-1][d-1]= coords(d, i);
			xrefgmx[i-1][d-1]=ref_coords(d,i);
		}
	}

	//compute rotation matrix
	calc_fit_R( natoms, weights_gmx, xrefgmx, xgmx, Rot );

	R.row_x( Vector( Rot[ 0 ][ 0 ], Rot[ 0 ][ 1 ], Rot[ 0 ][ 2 ] ) );
	R.row_y( Vector( Rot[ 1 ][ 0 ], Rot[ 1 ][ 1 ], Rot[ 1 ][ 2 ] ) );
	R.row_z( Vector( Rot[ 2 ][ 0 ], Rot[ 2 ][ 1 ], Rot[ 2 ][ 2 ] ) );

	for ( Size i = 1; i<=natoms ; i++ ) {
		Vector x( coords( 1, i ), coords( 2, i), coords( 3, i ) );
		x = R * x;
		coords( 1, i ) = x( 1 );
		coords( 2, i ) = x( 2 );
		coords( 3, i ) = x( 3 );
	}

	delete[] xgmx;
	delete[] xrefgmx;
	delete[] weights_gmx;
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); \
  a[k][l]=h+s*(g-h*tau);
#define DIM6 6
#define XX 0
#define YY 1
#define ZZ 2

void oprod(const rvec a,const rvec b,rvec c)
{
	c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
	c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
	c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}


void calc_fit_R(int natoms, Real const* weights, rvec const* xref, rvec const*x, matrix R)
{

	int    c,r,n,j,i,irot;
	double omega[ DIM6 ][ DIM6 ];
	double om[ DIM6 ] [ DIM6 ];
	double d[ DIM6 ],xnr,xpc;
	matrix vh,vk,u;
	Real   mn;
	int    index;
	Real   max_d;

	for ( i=0; i<DIM6; i++ ) {
		d[i]=0;
		for ( j=0; j<DIM6; j++ ) {
			omega[i][j]=0;
			om[i][j]=0;
		}
	}

	/* clear matrix U */
	for ( auto & i : u ) {
		for (double & j : i) j=0;
	}

	/*calculate the matrix U*/
	for ( n=0; n<natoms; n++ ) {
		if ( (mn = weights[ n ]) != 0.0 ) {
			for ( c=0; c<DIM; c++ ) {
				xpc=xref[n][c];
				for ( r=0; (r<DIM); r++ ) {
					xnr=x[n][r];
					u[c][r]+=mn*xnr*xpc;
				}
			}
		}
	}
	dump_matrix(DIM, u, tr);
	/*construct omega*/
	/*omega is symmetric -> omega==omega' */
	for ( r=0; r<DIM6; r++ ) {
		for ( c=0; c<=r; c++ ) {
			if ( r>=DIM && c<DIM ) {
				omega[r][c]=u[r-DIM][c];
				omega[c][r]=u[r-DIM][c];
			} else {
				omega[r][c]=0;
				omega[c][r]=0;
			}
		}
	}
	dump_matrix(DIM6, omega, tr);
	/*determine h and k*/
	jacobi( omega,d,om,&irot);
	/*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
	*int     natoms = number of rows and columns
	*real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
	*real       **v = v[0..n-1][0..n-1] contains the vectors in columns
	*int      *irot = number of jacobi rotations
	*/
	dump_matrix( 2*DIM, omega, tr );
	dump_matrix ( 2*DIM, om, tr );
	index=0; /* For the compiler only */

	/* Copy only the first two eigenvectors */
	for ( j=0; j<2; j++ ) {
		max_d=-1000;
		for ( i=0; i<DIM6; i++ ) {
			if ( d[i]>max_d ) {
				max_d=d[i];
				index=i;
			}
		}
		d[index]=-10000;
		for ( i=0; i<DIM; i++ ) {
			vh[j][i]=sqrt(2.0)*om[i][index];
			vk[j][i]=sqrt(2.0)*om[i+DIM][index];
		}
	}
	/* Calculate the last eigenvector as the outer-product of the first two.
	* This insures that the conformation is not mirrored and
	* prevents problems with completely flat reference structures.
	*/

	dump_matrix( DIM, vh, tr );
	dump_matrix( DIM, vk, tr );
	oprod(vh[0],vh[1],vh[2]);
	oprod(vk[0],vk[1],vk[2]);
	dump_matrix( DIM, vh, tr );
	dump_matrix( DIM, vk, tr );

	/*determine R*/
	for ( r=0; r<DIM; r++ ) {
		for ( c=0; c<DIM; c++ ) {
			R[r][c] = vk[0][r]*vh[0][c] +
				vk[1][r]*vh[1][c] +
				vk[2][r]*vh[2][c];
		}
	}
	dump_matrix( DIM, R, tr );
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); \
  a[k][l]=h+s*(g-h*tau);
#define DIM6 6
#define XX 0
#define YY 1
#define ZZ 2

void jacobi(double a[6][6],double d[],double v[6][6],int *nrot)
{
	int j,i;
	int iq,ip;
	double tresh,theta,tau,t,/*sm,*/s,h,g,c;
	double b[DIM6];
	double z[DIM6];
	int const n( DIM6 );
	for ( ip=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++ ) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for ( ip=0; ip<n; ip++ ) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for ( i=1; i<=50; i++ ) {
		double sm=0.0;
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				sm += fabs(a[ip][iq]);
			}
		}
		if ( sm == 0.0 ) {
			return;
		}
		if ( i < 4 ) {
			tresh=0.2*sm/(n*n);
		} else {
			tresh=0.0;
		}
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				g=100.0*fabs(a[ip][iq]);
				if ( i > 4 && fabs(d[ip])+g == fabs(d[ip])
						&& fabs(d[iq])+g == fabs(d[iq]) ) {
					a[ip][iq]=0.0;
				} else if ( fabs(a[ip][iq]) > tresh ) {
					h=d[iq]-d[ip];
					if ( fabs(h)+g == fabs(h) ) {
						t=(a[ip][iq])/h;
					} else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if ( theta < 0.0 ) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for ( j=0; j<ip; j++ ) {
						ROTATE(a,j,ip,j,iq);
					}
					for ( j=ip+1; j<iq; j++ ) {
						ROTATE(a,ip,j,j,iq);
					}
					for ( j=iq+1; j<n; j++ ) {
						ROTATE(a,ip,j,iq,j);
					}
					for ( j=0; j<n; j++ ) {
						ROTATE(v,j,ip,j,iq);
					}
					++(*nrot);
				}
			}
		}
		for ( ip=0; ip<n; ip++ ) {
			b[ip] +=  z[ip];
			d[ip]  =  b[ip];
			z[ip]  =  0.0;
		}
	}
	runtime_assert(0);
}


}
}
