// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief class to compute projection of a protein structure to principal component (PCA) eigenvectors ( as defined in file )
/// @author Oliver Lange


#if (defined _WIN32) && (!defined WIN_PYROSETTA)
#define ZLIB_WINAPI  // REQUIRED FOR WINDOWS
#endif

// Unit Headers
#include <protocols/evaluation/PCA.hh>

// Package Headers

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

//// C++ headers
#include <string>
#include <iostream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

/// @details Auto-generated virtual destructor
PCA::~PCA() = default;

static THREAD_LOCAL basic::Tracer tr( "protocols.evaluation.PCA", basic::t_info );

using namespace core;


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

void PCA::rotate_vec(int natoms,rvec *x,matrix R)
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

void PCA::add_vec( int natoms, rvec *x, rvec transvec ) {
	for ( int i=0; i<natoms; i++ ) {
		for ( int j = 0; j< DIM; j++ ) {
			x[i][j]+=transvec[j];
		}
	}
}


/// @brief read definition of PCA from file
void PCA::read_eigvec_file( std::string fn, pose::Pose const& pose,int nvec) {
	utility::io::izstream data( fn.c_str() );
	if ( !data ) {
		std::cerr << "ERROR:: Unable to open PCA file: "
			<< fn << std::endl;
		std::exit( 1 );
	}

	std::string line;
	getline(data,line); // header line
	getline(data,line); //nfit....
	std::istringstream line_stream ( line );
	std::string tag1, tag2, tag3;

	line_stream >> tag1 >> nfit_ >> tag2 >> npca_ >> tag3 >> nvec_;
	if ( nvec > 0 ) nvec_ = nvec;
	xref_.dimension( 3, nfit_ );
	xav_.dimension( 3, npca_ );
	ifit_.resize( nfit_ );
	ipca_.resize( npca_ );
	eigvec_.dimension( 3, npca_, nvec_ );

	getline(data,line); //AVERAGE
	if ( line != "AVERAGE" ) utility_exit_with_message(" tag AVERAGE missing ");
	read_structure( data, pose, ipca_, xav_,  "REFERENCE" );
	read_structure( data, pose, ifit_, xref_, "VECTORS" );
	for ( Size i=1; i<=nvec_; i++ ) {
		for ( Size k=1; k<=npca_; k++ ) {
			getline(data, line);
			std::istringstream line_stream( line );
			line_stream >> eigvec_( 1, k, i) >> eigvec_( 2, k, i) >> eigvec_(3, k, i);
		}
	}
}


/// @brief helper to read_eigvec_file: reads a protein structure from input file
void PCA::read_structure (
	std::istream& data,
	pose::Pose const& pose,
	IndexVector& ind,
	CoordVector& x,
	std::string endtag )
{
	std::string line;
	getline( data, line);
	int ct=1;
	while ( line != endtag ) {
		std::istringstream line_stream( line );
		std::string atomname;
		Size resnr;
		Real x1, x2, x3;
		line_stream >> atomname >> resnr >> x1 >> x2 >> x3;
		//---
		tr.Debug << "read PCA: " << atomname << " " << resnr << " " << x1 << " " << x2 << " " << x3 << "\n";
		ind[ ct ]=id::AtomID( pose.residue_type( resnr ).atom_index( atomname ), resnr );
		x( 1, ct ) = x1*10; x(2, ct) = x2*10; x(3, ct) = x3*10;
		//-
		getline( data, line);
		ct ++;
	}
}

void PCA::reset_x( Size n, CoordVector& x, rvec transvec ) {
	Size dim( 3 );
	// align center of mass to origin
	for ( Size k = 1; k <= dim; ++k ) {
		Real temp1 = 0.0;
		for ( Size j = 1; j <= n; ++j ) {
			temp1 += x(k,j);
		}
		temp1 /= 1.0*n;
		transvec[k-1]=-temp1;
		for ( Size j = 1; j <= n; ++j ) {
			x(k,j) -= temp1;
		}
	}
}

/// @brief compute projections for given pose
void PCA::eval( pose::Pose const& pose, ProjectionVector& proj ) {
	//fill Farray for fit
	rvec* xgmx;
	rvec* xrefgmx;
	rvec transvec;
	runtime_assert ( npca_ == nfit_ );// different fit- and analysis group doesn't work yet. some stupid bug.
	xgmx = new rvec[ npca_ ];
	xrefgmx = new rvec[ npca_ ];
	matrix Rot;
	CoordVector x;
	fill_coordinates( pose, ifit_, x );

	reset_x( nfit_, x, transvec );

	//transfer into C-style arrays
	for ( Size i = 1; i<=nfit_ ; i++ ) {
		for ( Size d = 1; d<=3; d++ ) {
			xgmx[i-1][d-1]= x(d, i)/10.0;
			xrefgmx[i-1][d-1]=xref_(d,i)/10.0;
		}
	}

	//compute rotation matrix
	calc_fit_R( nfit_, xrefgmx, xgmx, Rot );
	dump_matrix( 3, Rot, tr );

	fill_coordinates( pose, ipca_, x );
	//transfer into C-style arrays
	for ( Size i = 1; i<=npca_ ; i++ ) {
		for ( int d = 1; d<=3; d++ ) {
			xgmx[i-1][d-1]= x(d, i);
		}
	}
	add_vec( npca_, xgmx, transvec );
	rotate_vec( npca_, xgmx, Rot );

	tr.Trace << "rotated and translated\n";
	for ( Size i = 1; i<=npca_ ; i++ ) {
		for ( Size d = 1; d<=3; d++ ) {
			x(d, i) = xgmx[i-1][d-1];
			tr.Trace << x(d, i)/10.0 << " ";
		}
		tr.Trace << "\n";
	}
	//Compute projection
	proj.resize( nvec_ );
	for ( Size v = 1; v <= nvec_; v++ ) {
		proj[ v ]=0;
		for ( Size k = 1; k <= npca_; k++ ) {
			for ( Size d = 1; d <= 3; d++ ) {
				proj[ v ]+= (x( d, k)-xav_( d, k)) * eigvec_( d, k, v)/10.0;
			}
		}
	}

	delete[] xgmx;
	delete[] xrefgmx;
}

// @brief dump stuff on screen
void PCA::show( std::ostream& os ) {
	os << nfit_ << " " << npca_ << " " << nvec_ << std::endl;
	os << "AVERAGE" << std::endl;
	for ( Size i = 1; i <= npca_ ; i++ ) {
		os << ipca_[i] << " " << xav_(1, i ) << " " << xav_(2,i ) << " " << xav_(3,i ) << std::endl;
	}

	for ( Size i = 1; i <= nfit_ ; i++ ) {
		os << xref_( 1, i ) << " " << xref_( 2, i ) << " " << xref_( 3, i ) << std::endl;
	}

	for ( Size k = 1; k <= nvec_ ; k++ ) {
		for ( Size i = 1; i <= npca_ ; i++ ) {
			os << eigvec_( 1, i, k ) << " " << eigvec_( 2, i, k ) << " " << eigvec_( 3, i, k ) << std::endl;
		}
	}

}


/// @brief helper of eval: get the coordinates of interest from pose
void PCA::fill_coordinates(
	pose::Pose const& pose,
	IndexVector const& ind,
	CoordVector & x
)
{
	int natoms = 1;
	x.redimension( 3, ind.size() );
	 for ( auto const & it : ind ) {
		PointPosition vec( pose.xyz( it ) );
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			x(k+1,natoms) = vec[k];
		}
		++natoms;
	}
}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); \
  a[k][l]=h+s*(g-h*tau);
#define DIM6 6
#define XX 0
#define YY 1
#define ZZ 2

void PCA::oprod(const rvec a,const rvec b,rvec c)
{
	c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
	c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
	c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}


void PCA::calc_fit_R(int natoms,rvec *xp,rvec const* x,matrix R)
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
		for ( int j=0; j<DIM; j++ ) i[j]=0;
	}
	/*calculate the matrix U*/
	for ( n=0; (n<natoms); n++ ) {
		if ( (mn = 1.0) != 0.0 ) {
			for ( c=0; (c<DIM); c++ ) {
				xpc=xp[n][c];
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


void PCA::jacobi(double a[6][6],double d[],double v[6][6],int *nrot)
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

} //evaluation
} //protocols

