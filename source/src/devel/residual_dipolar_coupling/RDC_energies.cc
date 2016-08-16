// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Residual Dipolar Coupling scoring
/// @brief Residual Dipolar Coupling scoring
/// @details
///   Contains currently: LoopModeler
///
///
/// @author Vatsan Raman

//devel headers
#include <devel/residual_dipolar_coupling/RDC_main.hh>
#include <devel/residual_dipolar_coupling/RDC_energies.hh>

//core headers
#include <core/types.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

//utility headers
#include <utility/exit.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

//Objexx headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/Fmath.hh>

//C++ headers
#include <string>
#include <map>
#include <vector>

#include <utility/vector1.hh>


using namespace core;
using namespace ObjexxFCL;

namespace devel {
namespace residual_dipolar_coupling {

void eval_dipolar(
	core::pose::Pose const & pose,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines
)
{

	core::Size const number_of_rows( All_RDC_lines.size() );
	core::Size const ORDERSIZE = { 5 };
	FArray2D< core::Real > A( number_of_rows, ORDERSIZE );
	FArray1D< core::Real > b( number_of_rows );
	FArray1D< core::Real > x( number_of_rows );
	FArray2D< core::Real > vec( 3, 3 );
	core::Real Azz, eta;
	bool reject = false;

	std::cout << "Do nothing eval_dipolar !" << pose.total_residue() << " " << number_of_rows << std::endl;

	assemble_datamatrix( pose, All_RDC_lines, A, b);
	for ( core::Size i = 1; i <= number_of_rows; ++i ) {
		std::cout << "Matrices A & b  " << A(i,1) << " " << A(i,2) << " " << A(i,3) << " " << A(i,4) << " " << A(i,5) << " " << b(i) << std::endl;
	}

	calc_ordermatrix( number_of_rows, ORDERSIZE, A, b, x, reject );
	if ( reject ) std::cout << "SET SCORE VALUE, FIX THIS LATER " << std::endl;
	calc_orderparam( x, vec, Azz, eta );
	calc_dipscore( A, x, b, All_RDC_lines, ORDERSIZE, Azz );


}

////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
void assemble_datamatrix(
	core::pose::Pose const & pose,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines,
	FArray2D< core::Real > & A,
	FArray1D< core::Real > & b
)
{

	// numeric::xyzVector< Real > N,H;
	numeric::xyzVector< Real > umn;

	utility::vector1< devel::residual_dipolar_coupling::RDC >::const_iterator it;
	core::Size nrow( 0 );
	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {

		++nrow;
		umn = pose.residue(it->res1()).atom("N").xyz() - pose.residue(it->res1()).atom("H").xyz();

		core::Real umn_x = umn.x()/it->fixed_dist();
		core::Real umn_y = umn.y()/it->fixed_dist();
		core::Real umn_z = umn.z()/it->fixed_dist();

		//filling matrix A
		A( nrow, 1 ) = umn_y*umn_y - umn_x*umn_x;
		A( nrow, 2 ) = umn_z*umn_z - umn_x*umn_x;
		A( nrow, 3 ) = 2.0*umn_x*umn_y;
		A( nrow, 4 ) = 2.0*umn_x*umn_z;
		A( nrow, 5 ) = 2.0*umn_z*umn_y;

		//filling matrix b
		b( nrow ) = it->Reduced_Jdipolar();

		std::cout << "nrow " << nrow << std::endl;
		//  std::cout << "Matrix A " << A(nrow,1) << " " << A(nrow,2) << " " << A(nrow,3) << " " << A(nrow,4) << " " << A(nrow,5) << std::endl;
	}

}

////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
void calc_ordermatrix(
	core::Size const & nrow,
	core::Size const & ORDERSIZE,
	FArray2D< core::Real > & A,
	FArray1D< core::Real > & b,
	FArray1D< core::Real > & x,
	bool & reject
)
{

	std::cout << "Do nothing calc_ordermatrix !" << std::endl;

	core::Real factor = { 1e-6 }; // cutoff factor for singular values in svd

	FArray2D< core::Real > U( nrow, ORDERSIZE );
	FArray1D< core::Real > w( ORDERSIZE ); // singular values
	FArray2D< core::Real > v( ORDERSIZE, ORDERSIZE );

	core::Real wmin, wmax; // min and max values of w
	core::Size sing; // number of singular values in w
	core::Real Sxx;

	for ( core::Size i = 1; i <= nrow; ++i ) { // copy A
		U(i,1) = A(i,1); // copy into U
		U(i,2) = A(i,2);
		U(i,3) = A(i,3);
		U(i,4) = A(i,4);
		U(i,5) = A(i,5);
	}

	//***************** MAKE CHANGES IF NECESSARY ****************
	//block of code that adds zeroes if the number of rows are less than 5
	//not sure if we need that. The number of data lines are definitely expected
	//to be greater than 5.

	svdcmp( U, nrow, ORDERSIZE, w, v );

	wmax = 0.0;
	for ( core::Size j = 1; j <= ORDERSIZE; ++j ) {
		if ( w(j) > wmax ) wmax = w(j);
	}
	wmin = wmax * factor;
	sing = 0;
	for ( core::Size j = 1; j <= ORDERSIZE; ++j ) {
		if ( w(j) < wmin ) {
			w(j) = 0.0;
			++sing;
		}
	}
	if ( sing > core::Size( std::abs( int( nrow ) - int( ORDERSIZE ) ) ) ) {
		std::cout << "SVD yielded a matrix singular above expectation " <<
			"in get_ordermatrix" << std::endl;
	}

	// find solution for exact dipolar values

	svbksb( U, w, v, nrow, ORDERSIZE, b, x );

	// x components: (Syy,Szz,Sxy,Sxz,Syz)
	// check for acceptable values

	reject = false;
	Sxx = -x(1) - x(2);
	if ( Sxx < -0.5 || Sxx > 1.0 ) reject = true; // Sxx
	if ( x(1) < -0.5 || x(1) > 1.0 ) reject = true; // Syy
	if ( x(2) < -0.5 || x(2) > 1.0 ) reject = true; // Szz
	if ( x(3) < -0.75 || x(3) > 0.75 ) reject = true; // Sxy
	if ( x(4) < -0.75 || x(4) > 0.75 ) reject = true; // Sxz
	if ( x(5) < -0.75 || x(5) > 0.75 ) reject = true; // Syz

	if ( reject ) {
		std::cout << "order matrix not physically meaningful" << std::endl;
		//  try with errors on dipolar values? map error?
		//  score = 0.0;
		return;
	}


}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void svdcmp(
	FArray2D< core::Real > & a,
	core::Size const & m,
	core::Size const & n,
	FArray1D< core::Real > & w,
	FArray2D< core::Real > & v
)
{

	std::cout << "Do nothing svdcmp ! " << std::endl;

	//U    USES pythag
	core::Size i,its,j,jj,k,l,nm;
	FArray1D< core::Real > rv1( n );
	core::Real anorm, c, f, g, h, s, scale, x, y, z;
	g = 0.0;
	scale = 0.0;
	anorm = 0.0;
	l = 0;
	nm = 0;
	for ( i = 1; i <= n; ++i ) {
		l = i+1;
		rv1(i) = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( i <= m ) {
			for ( k = i; k <= m; ++k ) {
				scale += std::abs(a(k,i));
			}
			if ( scale != 0.0 ) {
				for ( k = i; k <= m; ++k ) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f = a(i,i);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,i) = f-g;
				for ( j = l; j <= n; ++j ) {
					s = 0.0;
					for ( k = i; k <= m; ++k ) {
						s += a(k,i)*a(k,j);
					}
					f = s/h;
					for ( k = i; k <= m; ++k ) {
						a(k,j) += f*a(k,i);
					}
				}
				for ( k = i; k <= m; ++k ) {
					a(k,i) *= scale;
				}
			}
		}
		w(i) = scale *g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( (i <= m) && (i != n) ) {
			for ( k = l; k <= n; ++k ) {
				scale += std::abs(a(i,k));
			}
			if ( scale != 0.0 ) {
				for ( k = l; k <= n; ++k ) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f = a(i,l);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,l) = f-g;
				for ( k = l; k <= n; ++k ) {
					rv1(k) = a(i,k)/h;
				}
				for ( j = l; j <= m; ++j ) {
					s = 0.0;
					for ( k = l; k <= n; ++k ) {
						s += a(j,k)*a(i,k);
					}
					for ( k = l; k <= n; ++k ) {
						a(j,k) += s*rv1(k);
					}
				}
				for ( k = l; k <= n; ++k ) {
					a(i,k) *= scale;
				}
			}
		}
		anorm = std::max(anorm,(std::abs(w(i))+std::abs(rv1(i))));
	}
	for ( i = n; i >= 1; --i ) {
		if ( i < n ) {
			if ( g != 0.0 ) {
				for ( j = l; j <= n; ++j ) {
					v(j,i) = (a(i,j)/a(i,l))/g;
				}
				for ( j = l; j <= n; ++j ) {
					s = 0.0;
					for ( k = l; k <= n; ++k ) {
						s += a(i,k)*v(k,j);
					}
					for ( k = l; k <= n; ++k ) {
						v(k,j) += s*v(k,i);
					}
				}
			}
			for ( j = l; j <= n; ++j ) {
				v(i,j) = 0.0;
				v(j,i) = 0.0;
			}
		}
		v(i,i) = 1.0;
		g = rv1(i);
		l = i;
	}
	for ( i = std::min(m,n); i >= 1; --i ) {
		l = i+1;
		g = w(i);
		for ( j = l; j <= n; ++j ) {
			a(i,j) = 0.0;
		}
		if ( g != 0.0 ) {
			g = 1.0/g;
			for ( j = l; j <= n; ++j ) {
				s = 0.0;
				for ( k = l; k <= m; ++k ) {
					s += a(k,i)*a(k,j);
				}
				f = (s/a(i,i))*g;
				for ( k = i; k <= m; ++k ) {
					a(k,j) += f*a(k,i);
				}
			}
			for ( j = i; j <= m; ++j ) {
				a(j,i) *= g;
			}
		} else {
			for ( j = i; j <= m; ++j ) {
				a(j,i) = 0.0;
			}
		}
		a(i,i) += 1.0;
	}
	for ( k = n; k >= 1; --k ) {
		for ( its = 1; its <= 30; ++its ) {
			for ( l = k; l >= 1; --l ) {
				nm = l-1;
				if ( (std::abs(rv1(l))+anorm) == anorm ) goto L2;
				if ( (std::abs(w(nm))+anorm) == anorm ) break;
			}
			c = 0.0;
			s = 1.0;
			for ( i = l; i <= k; ++i ) {
				f = s*rv1(i);
				rv1(i) *= c;
				if ( (std::abs(f)+anorm) == anorm ) break;
				g = w(i);
				h = pythag(f,g);
				w(i) = h;
				h = 1.0/h;
				c = (g*h);
				s = -(f*h);
				for ( j = 1; j <= m; ++j ) {
					y = a(j,nm);
					z = a(j,i);
					a(j,nm) = (y*c)+(z*s);
					a(j,i) = -(y*s)+(z*c);
				}
			}
			L2:
			z = w(k);
			if ( l == k ) {
				if ( z < 0.0 ) {
					w(k) = -z;
					for ( j = 1; j <= n; ++j ) {
						v(j,k) = -v(j,k);
					}
				}
				break;
			}
			if ( its == 30 ) utility_exit_with_message("no convergence in svdcmp \n" );
			x = w(l);
			nm = k-1;
			y = w(nm);
			g = rv1(nm);
			h = rv1(k);
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c = 1.0;
			s = 1.0;
			for ( j = l; j <= nm; ++j ) {
				i = j+1;
				g = rv1(i);
				y = w(i);
				h = s*g;
				g *= c;
				z = pythag(f,h);
				rv1(j) = z;
				c = f/z;
				s = h/z;
				f = (x*c)+(g*s);
				g = -(x*s)+(g*c);
				h = y*s;
				y *= c;
				for ( jj = 1; jj <= n; ++jj ) {
					x = v(jj,j);
					z = v(jj,i);
					v(jj,j) = (x*c)+(z*s);
					v(jj,i) = -(x*s)+(z*c);
				}
				z = pythag(f,h);
				w(j) = z;
				if ( z != 0.0 ) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = (c*g)+(s*y);
				x = -(s*g)+(c*y);
				for ( jj = 1; jj <= m; ++jj ) {
					y = a(jj,j);
					z = a(jj,i);
					a(jj,j) = (y*c)+(z*s);
					a(jj,i) = -(y*s)+(z*c);
				}
			}
			rv1(l) = 0.0;
			rv1(k) = f;
			w(k) = x;
		}
	}


}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
core::Real pythag(
	core::Real const & a,
	core::Real const & b
)
{

	core::Real pythag;

	core::Real absa = std::abs(a);
	core::Real absb = std::abs(b);
	if ( absa > absb ) {
		core::Real const ratio = absb/absa;
		pythag = absa * std::sqrt( 1.0 + ( ratio * ratio ) );
	} else {
		if ( absb == 0.0 ) {
			pythag = 0.0;
		} else {
			core::Real const ratio = absa/absb;
			pythag = absb * std::sqrt( 1.0 + ( ratio * ratio ) );
		}
	}
	return pythag;


}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void svbksb(
	FArray2D< core::Real > const & u,
	FArray1D< core::Real > const & w,
	FArray2D< core::Real > const & v,
	core::Size const & m,
	core::Size const & n,
	FArray1D< core::Real > const & b,
	FArray1D< core::Real > & x
)
{

	FArray1D< core::Real > tmp( n );
	core::Real s;

	for ( core::Size j = 1; j <= n; ++j ) {
		s = 0.0;
		if ( w(j) != 0.0 ) {
			for ( core::Size i = 1; i <= m; ++i ) {
				s += u(i,j) * b(i);
			}
			s /= w(j);
		}
		tmp(j) = s;
	}
	for ( core::Size j = 1; j <= n; ++j ) {
		s = 0.0;
		for ( core::Size jj = 1; jj <= n; ++jj ) {
			s += v(j,jj) * tmp(jj);
		}
		x(j) = s;
	}

}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void calc_orderparam(
	FArray1D< core::Real > x,
	FArray2D< core::Real > vec,
	core::Real & Azz,
	core::Real & eta
)
{

	using numeric::xyzMatrix;
	std::cout << "Do nothing calc_orderparam ! " << std::endl;

	FArray1D< core::Size > sort( 3 ); // sorted index to val, val(sort(1)) = largest abs val.
	core::Real temp1, temp2;

	// Assemble order matrix
	numeric::xyzMatrix< core::Real > S = numeric::xyzMatrix< core::Real >::rows( -x(1) - x(2), x(3), x(4),x(3), x(1), x(5), x(4), x(5), x(2) );


	numeric::xyzVector< core::Real > val; // Eigenvalues
	numeric::xyzMatrix< core::Real > xyz_vec; // Eigenvectors

	// Find eigenvalues, eigenvectors of symmetric matrix S
	val = eigenvector_jacobi( S, 1E-9, xyz_vec );

	// sort eigenvalues
	sort(1) = 1;
	sort(2) = 2;
	sort(3) = 3;

	if ( std::abs(val(1)) < std::abs(val(2)) ) {     // largest absolute value
		sort(2) = 1;
		sort(1) = 2;
	}
	if ( std::abs(val(sort(2))) < std::abs(val(3)) ) {
		sort(3) = sort(2);
		sort(2) = 3;
		if ( std::abs(val(sort(1))) < std::abs(val(3)) ) {
			sort(2) = sort(1);
			sort(1) = 3;
		}
	}


	Azz = val(sort(1));
	eta = (2.0/3.0) * std::abs(val(sort(2))-val(sort(3))/Azz);

	// sort eigen values      // largest to smallest : Azz,Ayy,Axx
	temp1 = val(sort(1));
	temp2 = val(sort(2));
	val(3) = val(sort(3));
	val(2) = temp2;
	val(1) = temp1;

	// sort eigen vectors
	for ( core::Size i = 1; i <= 3; ++i ) {
		temp1 = xyz_vec(i,sort(3));
		temp2 = xyz_vec(i,sort(2));
		vec(i,3) = xyz_vec(i,sort(1));
		vec(i,1) = temp1;
		vec(i,2) = temp2;
	}

	Azz = val(1);
	eta = (2.0/3.0) * (val(3)-val(2))/val(1);


}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////

void calc_dipscore(
	FArray2D< core::Real > const & A,
	FArray1D< core::Real > const & x,
	FArray1D< core::Real > const & b,
	utility::vector1<devel::residual_dipolar_coupling::RDC> const & All_RDC_lines,
	core::Size const & ORDERSIZE,
	core::Real const & Azz
)
{

	assert( Azz != 0 );

	core::Real score( 0.0 );

	utility::vector1< devel::residual_dipolar_coupling::RDC >::const_iterator it;
	core::Size nrow( 0 );
	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
		core::Real Jcalc( 0.0 );
		++nrow;
		for ( core::Size j = 1; j <= ORDERSIZE; ++j ) {
			Jcalc += A( nrow, j )*x(j);
		}
		score += ( b( nrow ) -Jcalc )*( b( nrow ) - Jcalc );

		std::cout << b(nrow)/it->invDcnst() << " " << Jcalc/it->invDcnst() << " " << numeric::square( b(nrow)/it->invDcnst() - Jcalc/it->invDcnst()) << " " << it->res1() << " " << it->res2() << std::endl;
	}

	std::cout << "score, nrow, Azz " << score << " " << nrow << " " << Azz << std::endl;

	score /= ( nrow*( Azz*Azz ) );

	std::cout << "Total score " << score << std::endl;

}

}//ResidualDipolarCoupling
}//devel
