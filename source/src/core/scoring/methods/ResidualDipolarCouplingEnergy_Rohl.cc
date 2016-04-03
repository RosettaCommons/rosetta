// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ResidualDipolarCouplingEnergy_Rohl.cc
/// @brief  RDC energy - comparing experimental RDC values to calculated values
/// @author Srivatsan Raman


//Unit headers
#include <core/scoring/methods/ResidualDipolarCouplingEnergy_Rohl.hh>
#include <core/scoring/methods/ResidualDipolarCouplingEnergy_RohlCreator.hh>
#include <core/scoring/ResidualDipolarCoupling_Rohl.hh>
#include <core/scoring/ResidualDipolarCoupling_Rohl.fwd.hh>

//Package headers

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

//utility headers
#include <utility/exit.hh>

//Objexx headers
#include <ObjexxFCL/Fmath.hh>

//C++ headers
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL;

/// @details This must return a fresh instance of the ResidualDipolarCouplingEnergy_Rohl class,
/// never an instance already in use
methods::EnergyMethodOP
ResidualDipolarCouplingEnergy_RohlCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ResidualDipolarCouplingEnergy_Rohl );
}

ScoreTypes
ResidualDipolarCouplingEnergy_RohlCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rdc_rohl );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCouplingEnergy_Rohl::ResidualDipolarCouplingEnergy_Rohl() :
	parent( EnergyMethodCreatorOP( new ResidualDipolarCouplingEnergy_RohlCreator ) )
{}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
EnergyMethodOP
ResidualDipolarCouplingEnergy_Rohl::clone() const
{
	return EnergyMethodOP( new ResidualDipolarCouplingEnergy_Rohl() );
}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const
{

	Real dipolar_score = eval_dipolar( pose );
	totals[ rdc_rohl ] = dipolar_score;

}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCoupling_Rohl const &
ResidualDipolarCouplingEnergy_Rohl::rdc_from_pose(
	pose::Pose & pose
) const
{
	ResidualDipolarCoupling_RohlOP rdc_info( retrieve_RDC_ROHL_from_pose( pose ) );
	if ( !rdc_info ) {
		rdc_info = ResidualDipolarCoupling_RohlOP( new ResidualDipolarCoupling_Rohl );
		store_RDC_ROHL_in_pose( rdc_info, pose );
	}
	return *rdc_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for RDC energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergy_Rohl::eval_dipolar(
	pose::Pose & pose
) const
{

	ResidualDipolarCoupling_Rohl const & rdc_data( rdc_from_pose( pose ) );
	utility::vector1< core::scoring::RDC_Rohl > All_RDC_lines( rdc_data.get_RDC_data() );

	Size const nrow( All_RDC_lines.size() ); //number of experimental couplins
	Size const ORDERSIZE = { 5 }; //Syy,Szz,Sxy,Sxz,Syz

	ObjexxFCL::FArray2D< Real > A( nrow, ORDERSIZE ); // N x 5, the 5 assymetric tensor elements per relevant vector in pose
	ObjexxFCL::FArray1D< Real > b( nrow );            // experimental values
	ObjexxFCL::FArray1D< Real > x( ORDERSIZE ); //previously dimensioned to nrow, which is wrong.
	ObjexxFCL::FArray1D< Real > weights( nrow ); //previously dimensioned to nrow, which is wrong.
	ObjexxFCL::FArray2D< Real > vec( 3, 3 );
	Real Azz, eta;
	bool reject = false;

	assemble_datamatrix( pose, All_RDC_lines, A, b, weights);
	/* for ( core::Size i = 1; i <= nrow; ++i ) {
	std::cout << "Matrices A & b  " << A(i,1) << " " << A(i,2) << " " << A(i,3) << " " << A(i,4) << " " << A(i,5) << " " << b(i) << std::endl;
	}
	*/
	calc_ordermatrix( nrow, ORDERSIZE, A, b, x, weights, reject );
	if ( reject ) std::cout << "SET SCORE VALUE, FIX THIS LATER " << std::endl;
	calc_orderparam( x, vec, Azz, eta );
	Real score( calc_dipscore( A, x, b, All_RDC_lines, ORDERSIZE, Azz )*nrow );

	return score;

}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::assemble_datamatrix(
	pose::Pose const & pose,
	utility::vector1< core::scoring::RDC_Rohl> const & All_RDC_lines,
	ObjexxFCL::FArray2D< Real > & A,
	ObjexxFCL::FArray1D< Real > & b,
	ObjexxFCL::FArray1D< Real > & weights
) const
{
	numeric::xyzVector< Real > umn;
	utility::vector1< core::scoring::RDC_Rohl >::const_iterator it;
	Size nrow( 0 );

	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
		++nrow;
		umn = pose.residue(it->res()).atom("N").xyz() - pose.residue(it->res()).atom("H").xyz();

		Real umn_x = umn.x()/it->fixed_dist();
		Real umn_y = umn.y()/it->fixed_dist();
		Real umn_z = umn.z()/it->fixed_dist();

		//filling matrix A
		A( nrow, 1 ) = umn_y*umn_y - umn_x*umn_x;
		A( nrow, 2 ) = umn_z*umn_z - umn_x*umn_x;
		A( nrow, 3 ) = 2.0*umn_x*umn_y;
		A( nrow, 4 ) = 2.0*umn_x*umn_z;
		A( nrow, 5 ) = 2.0*umn_z*umn_y;

		//filling matrix b
		b( nrow ) = it->Reduced_Jdipolar();
		weights( nrow ) = it->weight();
		//   std::cout << "nrow " << nrow << std::endl;
		//  std::cout << "Matrix A " << A(nrow,1) << " " << A(nrow,2) << " " << A(nrow,3) << " " << A(nrow,4) << " " << A(nrow,5) << std::endl;
	}
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::calc_ordermatrix(
	Size const & nrow,
	Size const & ORDERSIZE,
	ObjexxFCL::FArray2D< Real > & A,
	ObjexxFCL::FArray1D< Real > & b,
	ObjexxFCL::FArray1D< Real > & x,
	ObjexxFCL::FArray1D< Real > & weights,
	bool & reject
) const
{
	core::Real factor = { 1e-6 }; // cutoff factor for singular values in svd

	ObjexxFCL::FArray2D< Real > U( nrow, ORDERSIZE );
	ObjexxFCL::FArray1D< Real > w( ORDERSIZE ); // singular values
	ObjexxFCL::FArray2D< Real > v( ORDERSIZE, ORDERSIZE );
	ObjexxFCL::FArray1D< Real > bweighted( nrow ); // singular values
	Real wmin, wmax; // min and max values of w
	Size sing; // number of singular values in w
	Real Sxx;

	// why not U=A; ? should even be faster!
	Size ct_align = 0;
	for ( core::Size i = 1; i <= nrow; ++i ) { // copy A
		if ( weights( i ) <= 0.000001 ) continue;

		++ct_align;
		U( ct_align, 1 ) = A(i,1) * weights( i ); // copy into U
		U( ct_align, 2 ) = A(i,2) * weights( i );
		U( ct_align, 3 ) = A(i,3) * weights( i );
		U( ct_align, 4 ) = A(i,4) * weights( i );
		U( ct_align, 5 ) = A(i,5) * weights( i );
		bweighted( ct_align  ) = b( i ) * weights( i );
	}

	svdcmp( U, ct_align, ORDERSIZE, w, v );

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
	if ( (int)sing > std::abs( int( ct_align ) - int( ORDERSIZE ) ) ) {
		std::cout << "SVD yielded a matrix singular above expectation " <<
			"in get_ordermatrix" << std::endl;
	}

	// find solution for exact dipolar values
	svbksb( U, w, v, ct_align, ORDERSIZE, bweighted, x );

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
void ResidualDipolarCouplingEnergy_Rohl::svdcmp(
	ObjexxFCL::FArray2D< Real > & a,
	Size const & row,
	Size const & row_length,
	ObjexxFCL::FArray1D< Real > & w,
	ObjexxFCL::FArray2D< Real > & v
) const
{
	//U    USES pythag
	Size i,its,j,jj,k,l,nm;
	ObjexxFCL::FArray1D< core::Real > rv1( row_length );
	Real anorm, c, f, g, h, s, scale, x, y, z;
	g = 0.0;
	scale = 0.0;
	anorm = 0.0;
	l = 0;
	nm = 0;
	for ( i = 1; i <= row_length; ++i ) {
		l = i+1;
		rv1(i) = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( i <= row ) {
			for ( k = i; k <= row; ++k ) {
				scale += std::abs(a(k,i));
			}
			if ( scale != 0.0 ) {
				for ( k = i; k <= row; ++k ) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f = a(i,i);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,i) = f-g;
				for ( j = l; j <= row_length; ++j ) {
					s = 0.0;
					for ( k = i; k <= row; ++k ) {
						s += a(k,i)*a(k,j);
					}
					f = s/h;
					for ( k = i; k <= row; ++k ) {
						a(k,j) += f*a(k,i);
					}
				}
				for ( k = i; k <= row; ++k ) {
					a(k,i) *= scale;
				}
			}
		}
		w(i) = scale * g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( (i <= row) && (i != row_length) ) {
			for ( k = l; k <= row_length; ++k ) {
				scale += std::abs(a(i,k));
			}
			if ( scale != 0.0 ) {
				for ( k = l; k <= row_length; ++k ) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f = a(i,l);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,l) = f-g;
				for ( k = l; k <= row_length; ++k ) {
					rv1(k) = a(i,k)/h;
				}
				for ( j = l; j <= row; ++j ) {
					s = 0.0;
					for ( k = l; k <= row_length; ++k ) {
						s += a(j,k)*a(i,k);
					}
					for ( k = l; k <= row_length; ++k ) {
						a(j,k) += s*rv1(k);
					}
				}
				for ( k = l; k <= row_length; ++k ) {
					a(i,k) *= scale;
				}
			}
		}
		anorm = std::max(anorm,(std::abs(w(i))+std::abs(rv1(i))));
	}
	v(i,i) = 1.0;
	g = rv1(i);
	l = i;
	for ( i = row_length-1; i >= 1; --i ) {
		if ( g != 0.0 ) {
			for ( j = l; j <= row_length; ++j ) {
				v(j,i) = (a(i,j)/a(i,l))/g;
			}
			for ( j = l; j <= row_length; ++j ) {
				s = 0.0;
				for ( k = l; k <= row_length; ++k ) {
					s += a(i,k)*v(k,j);
				}
				for ( k = l; k <= row_length; ++k ) {
					v(k,j) += s*v(k,i);
				}
			}
		}
		for ( j = l; j <= row_length; ++j ) {
			v(i,j) = 0.0;
			v(j,i) = 0.0;
		}
		v(i,i) = 1.0;
		g = rv1(i);
		l = i;
	}

	for ( i = std::min(row, row_length); i >= 1; --i ) {
		l = i+1;
		g = w(i);
		for ( j = l; j <= row_length; ++j ) {
			a(i,j) = 0.0;
		}
		if ( g != 0.0 ) {
			g = 1.0/g;
			for ( j = l; j <= row_length; ++j ) {
				s = 0.0;
				for ( k = l; k <= row; ++k ) {
					s += a(k,i)*a(k,j);
				}
				f = (s/a(i,i))*g;
				for ( k = i; k <= row; ++k ) {
					a(k,j) += f*a(k,i);
				}
			}
			for ( j = i; j <= row; ++j ) {
				a(j,i) *= g;
			}
		} else {
			for ( j = i; j <= row; ++j ) {
				a(j,i) = 0.0;
			}
		}
		a(i,i) += 1.0;
	}

	for ( k = row_length; k >= 1; --k ) {
		for ( its = 1; its <= 30; ++its ) {
			bool skipnow(false);
			for ( l = k; l >= 1; --l ) {
				nm = l-1;
				if ( (std::abs(rv1(l))+anorm) == anorm ) {
					skipnow=true;
					break;
				}
				if ( (std::abs(w(nm))+anorm) == anorm ) break;
			}
			if ( !skipnow ) {
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
					for ( j = 1; j <= row; ++j ) {
						y = a(j,nm);
						z = a(j,i);
						a(j,nm) = (y*c)+(z*s);
						a(j,i) = -(y*s)+(z*c);
					}
				}
			} //if(!skipnow)
			z = w(k);
			if ( l == k ) {
				if ( z < 0.0 ) {
					w(k) = -z;
					for ( j = 1; j <= row_length; ++j ) {
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
				for ( jj = 1; jj <= row_length; ++jj ) {
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
				for ( jj = 1; jj <= row; ++jj ) {
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
Real ResidualDipolarCouplingEnergy_Rohl::pythag(
	Real const & a,
	Real const & b
) const
{
	Real pythag;

	Real absa = std::abs(a);
	Real absb = std::abs(b);
	if ( absa > absb ) {
		Real const ratio = absb/absa;
		pythag = absa * std::sqrt( 1.0 + ( ratio * ratio ) );
	} else {
		if ( absb == 0.0 ) {
			pythag = 0.0;
		} else {
			Real const ratio = absa/absb;
			pythag = absb * std::sqrt( 1.0 + ( ratio * ratio ) );
		}
	}
	return pythag;
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::svbksb(
	ObjexxFCL::FArray2D< Real > const & u,
	ObjexxFCL::FArray1D< Real > const & w,
	ObjexxFCL::FArray2D< Real > const & v,
	Size const & m,
	Size const & n,
	ObjexxFCL::FArray1D< Real > const & b,
	ObjexxFCL::FArray1D< Real > & x
) const
{
	ObjexxFCL::FArray1D< core::Real > tmp( n );
	Real s;

	for ( Size j = 1; j <= n; ++j ) {
		s = 0.0;
		if ( w(j) != 0.0 ) {
			for ( Size i = 1; i <= m; ++i ) {
				s += u(i,j) * b(i);
			}
			s /= w(j);
		}
		tmp(j) = s;
	}
	for ( Size j = 1; j <= n; ++j ) {
		s = 0.0;
		for ( Size jj = 1; jj <= n; ++jj ) {
			s += v(j,jj) * tmp(jj);
		}
		x(j) = s;
	}
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::calc_orderparam(
	ObjexxFCL::FArray1D< Real > x,
	ObjexxFCL::FArray2D< Real > vec,
	Real & Azz,
	Real & eta
) const
{
	using numeric::xyzMatrix;

	ObjexxFCL::FArray1D< Size > sort( 3 ); // sorted index to val, val(sort(1)) = largest abs val.
	Real temp1, temp2;

	// Assemble order matrix
	numeric::xyzMatrix< Real > S = numeric::xyzMatrix< core::Real >::rows( -x(1) - x(2), x(3), x(4),x(3), x(1), x(5), x(4), x(5), x(2) );


	numeric::xyzVector< Real > val; // Eigenvalues
	numeric::xyzMatrix< Real > xyz_vec; // Eigenvectors

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
	for ( Size i = 1; i <= 3; ++i ) {
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
Real ResidualDipolarCouplingEnergy_Rohl::calc_dipscore(
	ObjexxFCL::FArray2D< Real > const & A,
	ObjexxFCL::FArray1D< Real > const & x,
	ObjexxFCL::FArray1D< Real > const & b,
	utility::vector1< core::scoring::RDC_Rohl > const & All_RDC_lines,
	Size const & ORDERSIZE,
	Real const & Azz
) const
{

	debug_assert( Azz != 0 );

	Real score( 0.0 );

	utility::vector1< core::scoring::RDC_Rohl >::const_iterator it;
	Size nrow( 0 );
	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
		Real Jcalc( 0.0 );
		++nrow;
		for ( Size j = 1; j <= ORDERSIZE; ++j ) {
			Jcalc += A( nrow, j )*x(j);
		}
		score += ( b( nrow ) -Jcalc )*( b( nrow ) - Jcalc ); //these are reduced Jd

		//v  std::cout << b(nrow)/it->invDcnst() << " " << Jcalc/it->invDcnst() << " " << numeric::square( b(nrow)/it->invDcnst() - Jcalc/it->invDcnst() ) << " " << it->res() << std::endl;
	}

	//std::cout << "score, nrow, Azz " << score << " " << nrow << " " << Azz << std::endl;

	score /= ( nrow*( Azz*Azz ) );

	//std::cout << "Total score " << score << std::endl;

	return score;

}
core::Size
ResidualDipolarCouplingEnergy_Rohl::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core
