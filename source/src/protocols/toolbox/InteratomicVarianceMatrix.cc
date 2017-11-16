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


#include <protocols/toolbox/InteratomicVarianceMatrix.hh>

#include <core/types.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1P.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>

//// C++ headers
#include <string>
#include <iostream>

#include <numeric/model_quality/RmsData.hh>
#include <ObjexxFCL/FArray3.hh>


namespace protocols {
namespace toolbox {

using namespace ObjexxFCL;

static basic::Tracer tr( "protocols.evaluation.PCA", basic::t_info );

using namespace core;
using namespace numeric::model_quality; //for rms functions

/// @brief
//compute matrix of distance variances:
// in: coords : 3 x natoms x ndecoys
// out: ivm_ : natoms x natoms with ivm(i,j) = VAR_n ( | x_i( n ) - x_j ( n ) | )  n==decoys, i,j=atoms
void InteratomicVarianceMatrix::init( Size atoms_in, Size n_decoys, ObjexxFCL::FArray3_double const& coords ) {
	n_atoms_ = atoms_in;
	ivm_.redimension( n_atoms_, n_atoms_, 0.0 );
	Real const invn( 1.0/n_decoys );
	tr.Debug << "IVM-Matrix \n";
	for ( Size i = 1; i <= n_atoms(); i++ ) {
		for ( Size j = i+1; j<=n_atoms(); j++ ) {
			Real var( 0.0 );
			Real av( 0.0 );

			for ( Size n = 1; n<=n_decoys; n++ ) {
				Vector xi( coords( 1, i, n ), coords( 2, i, n ), coords( 3, i, n ));
				Vector xj( coords( 1, j, n ), coords( 2, j, n ), coords( 3, j, n ));
				Real dist = xi.distance(xj);
				var += dist * dist * invn;
				av += dist * invn;
			}
			ivm_( j, i ) = var - av*av;
			ivm_( i, j ) = ivm_( j, i);
			tr.Debug << i << " " << j << " " << ivm_( i, j ) << "\n";
		}
	}
	tr.Debug << std::endl;

	utility::io::ozstream os("ivm.dat");
	for ( Size i = 1; i <= n_atoms(); i++ ) {
		for ( Size j = 1; j<=n_atoms(); j++ ) os << ivm_( i, j ) << " ";
		os << std::endl;
	}
}


/// @brief compute order parameter for atom i, defined as number of j atoms whose ivm(i,j)<epsilon^2
void InteratomicVarianceMatrix::order_parameter( Real epsilon, ObjexxFCL::FArray1_double& T ) {
	Real const epsilon2( epsilon*epsilon );
	tr.Debug << "T( " << epsilon <<  " )\n";
	for ( Size i = 1; i<=n_atoms(); i++ ) {
		T( i ) = 0.0;
		Real invn = 1.0/n_atoms();
		for ( Size j = 1; j<=n_atoms(); j++ ) {
			if ( ivm_( i, j) < epsilon2 ) T( i ) += invn;
		}
		tr.Debug << i << " " << T( i ) << "\n";
	}
	tr.Debug << std::endl;
}

/// @brief compute order parameter for atom i, defined as number of j atoms whose ivm(i,j)<epsilon^2
Real InteratomicVarianceMatrix::kurtosis( ObjexxFCL::FArray1_double& T ) {
	Real invn = 1.0/n_atoms();
	Real second_moment = 0.0;
	Real fourth_moment = 0.0;
	Real av = 0.0;
	for ( Size j = 1; j<=n_atoms(); j++ ) {
		av += T( j )*invn;
	}
	for ( Size j = 1; j<=n_atoms(); j++ ) {
		Real x = T( j ) - av;
		Real x2 = x*x;
		second_moment += x2*invn;
		fourth_moment += x2 * x2 *invn;
	}
	return  fourth_moment / (second_moment*second_moment);
}

void InteratomicVarianceMatrix::optimize_kurtosis( Size ngrid, Real lb, Real ub ) {
	Real depsilon =  ( ub - lb ) / ngrid;
	Real epsilon = lb;
	FArray1D_double grid( ngrid, 0.0 );
	FArray1D_double kurt( ngrid, 0.0 );

	FArray2D_double T( n_atoms(), ngrid, 0.0 );
	tr.Info << "kurtosis computed:\n";
	for ( Size i = 1; i <= ngrid; i++ ) {
		FArray1P_double Tslice( T( 1, i ), n_atoms() );
		grid( i ) = epsilon;
		order_parameter( epsilon, Tslice );
		kurt( i ) = kurtosis( Tslice );
		tr.Info << epsilon << " " << kurt( i ) << "\n";
		epsilon += depsilon;
	}
	tr.Info << std::endl;
	{
		utility::io::ozstream os("order.dat");
		for ( Size j = 1; j<= n_atoms(); j++ ) {
			for ( Size i = 1; i<= ngrid; i++ )  os << T( j, i ) << " ";
			os << "\n";
		}
	}
	{
		utility::io::ozstream os("kurt.dat");
		for ( Size i = 1; i<= ngrid; i++ )  os << grid( i ) << " " << kurt( i ) << "\n";
	}
}

} //evaluation
} //protocols


