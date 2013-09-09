// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/normalmode/NormalMode.cc
/// @brief   
/// @detailed
/// @author  Hahnbeom Park

// Unit headers

// Package headers
#include <protocols/normalmode/NormalMode.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix_operations.hh>
//#include <numeric/xyzMatrix.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/dssp/Dssp.hh>

/// C++ headers
#include <string>
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

static basic::Tracer TR("protocols.normalmode.NormalMode");

namespace protocols {
namespace normalmode {

using namespace core;
using namespace protocols::normalmode;

NormalMode::NormalMode()
{
	use_uniform_k_ = true;
	k_uniform_ = 1.0;
	torsion_ = false;
	mode_ = "CA";
	dist2_cut_ = 100.0;
	use_phi_ = false;
	use_psi_ = false;
	eckart_correction_ = false;
	torsions_using_assigned_ = false;
}

NormalMode::NormalMode( std::string const mode,
												Real const distcut )
{
	use_uniform_k_ = true;
	k_uniform_ = 1.0;
	torsion_ = false;
	mode_ = mode;
	dist2_cut_ = distcut*distcut;
	use_phi_ = false;
	use_psi_ = false;
	eckart_correction_ = false;
	torsions_using_assigned_ = false;
}

NormalMode::~NormalMode(){}

void
NormalMode::torsion( bool const torsion, 
										 bool const use_phi, 
										 bool const use_psi,
										 bool const eckart_correction ){ 
	torsion_ = torsion;
	use_phi_ = use_phi;
	use_psi_ = use_psi;
	eckart_correction_ = eckart_correction;
}

void
NormalMode::set_harmonic_constants( Real const & k_uniform ){
	k_uniform_ = k_uniform;
}

void
NormalMode::set_harmonic_constants( Real const & k_short,
																		Real const & k_SS,
																		Real const & k_long ){
	use_uniform_k_ = false;
	k_short_ = k_short;
	k_SS_ = k_SS;
	k_long_ = k_long;

}

void 
NormalMode::solve( pose::Pose const &pose )
{	
	if( torsion() ){
		TR << "Solving Torsional Normal Mode (phi/psi:" << use_phi_ << "/" << use_psi_ << "), ";
	} else {
		TR << "Solving Cartesian Noraml Mode, ";
	}
	TR << "Distance cut set as " << std::sqrt(dist2_cut_);
	if( !use_uniform_k_ ){
		TR << ", using 3-state spring constants:" << k_short_ << " " << k_SS_ << " " << k_long_ << std::endl;
	} else {
		TR << ", using uniform spring constants:" << k_uniform_ << std::endl;
	}

	prepare_coord( pose );

	set_harmonic_constant_map( pose );

	utility::vector1< utility::vector1< Real > > U, U_cart;

	Size n;
	if( torsion() ){
		// New way
		U = make_Hessian_TNM();

		// Old way
		//U_cart = make_Hessian_ANM();
		//U = convert_to_torsion_crd( U_cart );

		n = ntor();
		eigvec_tor_.resize( n );

	} else {
		U = make_Hessian_ANM();
		n = 3*natm();
		eigvec_cart_.resize( n );
	}

	// Allocate size given degree of freedom
	eigval_.resize( n );
	utility::vector1< utility::vector1< Real > > eigvec_tmp( n );
	for( Size i = 1; i <= n; ++i ) eigvec_tmp[i].resize( n );

	// Solve eigenvalue problem
	svdcmp( U, n, n, eigval_, eigvec_tmp );
	TR << "Diagonalization done." << std::endl;

	// Sort matrix in descending order of eigenvalues
	eigsrt( U, eigval_ );

	if( torsion() ){
		eigvec_tor_ = U;

	} else { 	// Cartesian only: Convert list into vector type
		for( Size i = 1; i <= 3*natm(); ++i ) eigvec_cart_[i].resize( natm() );

		for( Size i=1; i<= 3*natm(); ++i ){
			for( Size j=1; j<= n; ++j ){ // 
				Size iatm = (j-1)/3+1;
				Size k = (j-1)%3;
				eigvec_cart_[i][iatm][k] = U[j][i];
			}
		}
	}

	// Convert eigen value info into importance
	importance_.resize( nmode() );
	Real inv_sum( 0.0 );
	for( Size i_mode = 1; i_mode <= nmode(); ++i_mode ){	inv_sum += 1.0/get_eigval( i_mode ); }

	for( Size i_mode = 1; i_mode <= nmode(); ++i_mode ) {
		importance_[i_mode] = (1.0/get_eigval(i_mode))/inv_sum;
	}

}

//-------------------------------------------------------------
// privates
//-------------------------------------------------------------

void
NormalMode::prepare_coord( pose::Pose const & pose ){

	atomID_.resize( 0 );
	xyz_.resize( 0 );
	s_.resize( 0 );
	e_.resize( 0 );
	tau_.resize( 0 );
	a_to_i_.resize( 0 );
	i_to_a_.resize( 0 );

	// Default - assign full residue list when the list is not specified
	if( !torsions_using_assigned_ ){
		for( Size ires = 1; ires <= pose.total_residue(); ++ires )
			torsions_using_.push_back( ires );
	}
	
	Vector xyzCA, xyzN, xyzC;
	Size n( 0 );

	//std::cout << "prepare_coord, nmode " << torsions_using_.size() << std::endl;

	for( Size i = 1; i <= torsions_using_.size(); ++i ){
		Size const ires = torsions_using_[i];
		//std::cout << "i/ires " << i << " " << ires << std::endl;
		conformation::Residue const &rsd( pose.residue( ires ) );

		for( Size iatm = 1; iatm <= rsd.natoms(); ++iatm ){
			std::string const atmname( rsd.atom_name( iatm ) );
			Vector crd = rsd.xyz( iatm );

			// Let's make only for CA for now
			if ( mode().compare( "CA" ) == 0 && atmname.compare( " CA " ) == 0 ){
				id::AtomID atmID( iatm, ires );
				atomID_.push_back( atmID );
				xyz_.push_back( crd );
				i_to_a_.push_back( n );
			}

			// Add atomID depending on representation
			if( atmname.compare( " CA " ) == 0 ){
				xyzCA = crd;

			} else if(atmname.compare( " N  " ) == 0 ){
				xyzN = crd;

			}	else if(atmname.compare( " C  " ) == 0 ){
				xyzC = crd;

			}
		}

		// Skip C-term residue since it does not give any effect to Calpha atoms
		if( rsd.is_upper_terminus() ) continue; 

		// Setup torsion index
		if( !rsd.is_lower_terminus() && use_phi_ ){ // not N-term
			s_.push_back( xyzN ); // origin for phi
			Vector dxyz( xyzCA - xyzN );
			dxyz /= std::sqrt( dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2] );
			e_.push_back( dxyz ); // vertex for phi
			tau_.push_back( -dxyz.cross( xyzN ) ); // shift vector

			a_to_i_.push_back( xyz_.size() );
			n = s_.size();
			core::id::TorsionID torID( rsd.seqpos(), core::id::BB, 1 );
			torID_.push_back( torID );

			TR.Debug << "PHI - a,ires,s,e,a_to_i: " << n << " " << ires;
			TR.Debug << " " << std::setw(8) << s_[n][0];
			TR.Debug << " " << std::setw(8) << s_[n][1];
			TR.Debug << " " << std::setw(8) << s_[n][2];
			TR.Debug << " " << std::setw(8) << e_[n][0];
			TR.Debug << " " << std::setw(8) << e_[n][1];
			TR.Debug << " " << std::setw(8) << e_[n][2];
			TR.Debug << " " << std::setw(8) << a_to_i_[n];
			TR.Debug << std::endl;
		}

		// not C-term nor N-term: N-term psi doesn't have downstream Calpha 
		// so Hessian will be always zero
		if( !rsd.is_upper_terminus() && !rsd.is_lower_terminus()
				&& use_psi_ ){ 
			s_.push_back( xyzCA ); // origin for psi
			Vector dxyz( xyzC - xyzCA );
			dxyz /= std::sqrt( dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2] );
			e_.push_back( dxyz ); // vertex for psi
			tau_.push_back( -dxyz.cross( xyzCA ) ); // shift vector

			a_to_i_.push_back( xyz_.size() );
			n = s_.size();
			core::id::TorsionID torID( rsd.seqpos(), id::BB, 2 );
			torID_.push_back( torID );

			TR.Debug << "PSI - a,ires,s,e,a_to_i: " << n << " " << ires;
			TR.Debug << " " << std::setw(8) << s_[n][0];
			TR.Debug << " " << std::setw(8) << s_[n][1];
			TR.Debug << " " << std::setw(8) << s_[n][2];
			TR.Debug << " " << std::setw(8) << e_[n][0];
			TR.Debug << " " << std::setw(8) << e_[n][1];
			TR.Debug << " " << std::setw(8) << e_[n][2];
			TR.Debug << " " << std::setw(8) << a_to_i_[n];
			TR.Debug << std::endl;
		}

	}

	// Correction for a_to_i: include C-term atom
	//a_to_i_[ a_to_i_.size() ] = xyz_.size();i
	
}

void
NormalMode::set_harmonic_constant_map( pose::Pose const &pose ){

	// Initialize
	k_.resize( natm() );
	for( Size iatm = 1; iatm <= natm(); ++iatm ){
		k_[iatm].resize( natm() );
	}

	// Get copy of pose for SS assignment
	pose::Pose pose_tmp( pose );
	core::scoring::dssp::Dssp dssp( pose_tmp );
	dssp.insert_ss_into_pose( pose_tmp );

	// Store consecutive secondary structure number
	utility::vector1< Size > SSno( pose.total_residue(), 0 );
	std::map< Size, Size > Hbond_map;
	Real const Hbond_d2( 3.5*3.5 );

	if( !use_uniform_k_ ){

		Size i_ss( 0 );
		char prvSS = 'X';
		for( Size ires = 1; ires <= pose_tmp.total_residue(); ++ires ){
			char const &SS_i = pose_tmp.secstruct( ires );

			Vector const &Ncrd = pose_tmp.residue( ires ).xyz( " N  " );

			//std::cout << "ires/SS: " << ires << " " <<  SS_i << std::endl;
			if( SS_i == 'L' ) continue; // Skip coil in any case

			// For consecutive 2ndary structure segment
			if( SS_i != prvSS ){ // If not the same SS as prv
				i_ss ++;
			}
			SSno[ires] = i_ss;
			prvSS = SS_i;

			// For Hbond pair
			Hbond_map[ires] = 0; // placeholder
			Real mind2( 100.0 );
			for( Size jres = 1; jres <= pose_tmp.total_residue(); ++jres ){
				char const &SS_j = pose_tmp.secstruct( jres );

				if( std::abs(int(ires - jres)) < 3. || SS_j == 'L' ) continue;

				Vector const &Ocrd = pose_tmp.residue( jres ).xyz( " O  " );
				Vector const vNO = Ncrd - Ocrd;
				Real const d2 = vNO.dot( vNO );

				if( d2 < Hbond_d2 && d2 < mind2 ){
					Hbond_map[ires] = jres;
					mind2 = d2;
				}
			}

		}
	}

	// Setup harmonic constants
	// Note that iatm, jatm represents C-alpha atoms 
	for( Size iatm = 1; iatm <= natm()-1; ++iatm ){
		Size ires = atomID_[iatm].rsd();
		for( Size jatm = iatm+1; jatm <= natm(); ++jatm ){
			Size jres = atomID_[jatm].rsd();

			if( use_uniform_k_ ){
				k_[iatm][jatm] = k_uniform_;
				k_[jatm][iatm] = k_uniform_;

			} else {
				if( abs(int(ires - jres)) <= 1 ){
					k_[iatm][jatm] = k_short_;
					k_[jatm][iatm] = k_short_;

					/* 
				 //version for consecutive 2ndary structure segment
				} else if ( SSno[ires] > 0 && SSno[jres] > 0 && SSno[ires] == SSno[jres] ){
					k_[iatm][jatm] = k_SS_;
					k_[jatm][iatm] = k_SS_;

					*/
				} else if( Hbond_map[iatm] == jatm || Hbond_map[jatm] == iatm ) {
					k_[iatm][jatm] = k_SS_;
					k_[jatm][iatm] = k_SS_;

				// version for residue pair connected by 
				} else {
					k_[iatm][jatm] = k_long_;
					k_[jatm][iatm] = k_long_;

				}

			}
		}
	}

}

utility::vector1< utility::vector1< Real > >
NormalMode::make_Hessian_ANM( ){
	
	// allocate array
	utility::vector1< utility::vector1< Real > > H( 3*natm() );
	for( Size i = 1; i <= 3*natm(); ++i ) H[i].resize( 3*natm() );

	// First, fill in off-diagonal elements
	for( Size i = 1; i <= natm()-1; ++i ){
		for( Size j = i+1; j <= natm(); ++j ){

			Vector const dxyz( xyz(i) - xyz(j) );
			Real dist2( dxyz.dot( dxyz ) );

			if( dist2 < dist2_cut() ){

				Real k =  get_k( i, j )/dist2;
				Real const dxy( dxyz[0]*dxyz[1] );
				Real const dyz( dxyz[1]*dxyz[2] );
				Real const dxz( dxyz[2]*dxyz[0] );
				H[3*j-2][3*i-2] = k*dxyz[0]*dxyz[0];
				H[3*j-2][3*i-1] = k*dxy;
				H[3*j-2][3*i  ] = k*dxz;
				H[3*j-1][3*i-2] = k*dxy;
				H[3*j-1][3*i-1] = k*dxyz[1]*dxyz[1];
				H[3*j-1][3*i  ] = k*dyz;
				H[3*j  ][3*i-2] = k*dxz;
				H[3*j  ][3*i-1] = k*dyz;
				H[3*j  ][3*i  ] = k*dxyz[2]*dxyz[2];

				// Copy it to j,i-th position
				H[3*i-2][3*j-2] = H[3*j-2][3*i-2];
				H[3*i-2][3*j-1] = H[3*j-2][3*i-1];
				H[3*i-2][3*j  ] = H[3*j-2][3*i  ];
				H[3*i-1][3*j-2] = H[3*j-1][3*i-2];
				H[3*i-1][3*j-1] = H[3*j-1][3*i-1];
				H[3*i-1][3*j  ] = H[3*j-1][3*i  ];
				H[3*i  ][3*j-2] = H[3*j  ][3*i-2];
				H[3*i  ][3*j-1] = H[3*j  ][3*i-1];
				H[3*i  ][3*j  ] = H[3*j  ][3*i  ];
			}
		}
	}

	// Get diagonal elements
	for( Size i = 1; i <= natm(); ++i ){
		for( Size j = 1; j <= natm(); ++j ){
			if( j != i ){
				H[3*i-2][3*i-2] -= H[3*i-2][3*j-2];
				H[3*i-2][3*i-1] -= H[3*i-2][3*j-1];
				H[3*i-2][3*i  ] -= H[3*i-2][3*j  ];
				H[3*i-1][3*i-2] -= H[3*i-1][3*j-2];
				H[3*i-1][3*i-1] -= H[3*i-1][3*j-1];
				H[3*i-1][3*i  ] -= H[3*i-1][3*j  ];
				H[3*i  ][3*i-2] -= H[3*i  ][3*j-2];
				H[3*i  ][3*i-1] -= H[3*i  ][3*j-1];
				H[3*i  ][3*i  ] -= H[3*i  ][3*j  ];
			}
		}
	}

	TR.Debug << "Hessian" << std::endl;
	for( Size i = 1; i <= 3*natm(); ++i ){
		for( Size j = 1; j <= 3*natm(); ++j ){
			TR.Debug << " " << std::setw(9) << H[i][j];
		}
		TR.Debug << std::endl;
	}

	return H;
}

utility::vector1< utility::vector1< Real > >
NormalMode::make_Hessian_TNM( ){
	
	// allocate array
	utility::vector1< utility::vector1< Real > > H( ntor() );
	for( Size i = 1; i <= ntor(); ++i ) H[i].resize( ntor() );

	// Correct Tau
	Vector Rcom( 0.0, 0.0, 0.0 );
	for( Size iatm = 1; iatm <= natm(); ++iatm ) Rcom += xyz( iatm );
	Rcom /= natm();

	for( Size k = 1; k <= ntor(); ++k )	tau_[k] += e_[k].cross( Rcom );

	utility::vector1< ContactStruct > cos;
	// First, get contact list
	for( Size i = 1; i <= natm()-1; ++i ){
		for( Size j = i+1; j <= natm(); ++j ){
			Vector const dxyz( xyz(i) - xyz(j) );
			Real dist2( dxyz.dot( dxyz ) );

			if( dist2 < dist2_cut() ){
				Real k =  get_k( i, j )/dist2;

				ContactStruct co;
				co.k = k;
				co.r21 = -dxyz; // r2 - r1
				co.c21 = (xyz(j) - Rcom).cross( xyz(i) - Rcom ); // r2 x r1
				co.d2 = dist2;
				co.tor1 = i_to_a( i );
				co.tor2 = i_to_a( j );
				cos.push_back( co );
			}
		}
	}

	// below is original code, but we can simplify?
	// Second, build eta & tau matrix
	utility::vector1< Vector > zero( ntor(), Vector( 0.0, 0.0, 0.0 ) );

	// Zero matrix
	utility::vector1< utility::vector1 < Vector > > eta_ia( ntor(), zero );
	utility::vector1< utility::vector1 < Vector > > tau_ia( ntor(), zero );

	// Fill upper diagonal
	for( Size k = 1; k <= ntor(); ++k ){
		for( Size l = k; l <= ntor(); ++l ){
			eta_ia[l][k] = e_[k];
			tau_ia[l][k] = tau_[k];
		}
	}

	// Then, fill in  following the rule:
	// given torsion a, b
	utility::vector1< Real > dk( ntor(), 0.0 );

	for( Size i_co = 1; i_co <= cos.size(); ++i_co ){
		ContactStruct const &co = cos[i_co];

		// Determine set of rotable axes
		// TODO: We need to make this generic for chain breaks 
		utility::vector1< Vector > tau1;
		utility::vector1< Vector > tau2;
		utility::vector1< Vector > eta1;
		utility::vector1< Vector > eta2;
		if( co.tor1 == 0 ){
			tau1 = zero;
			eta1 = zero;
		} else {
			tau1 = tau_ia[co.tor1];
			eta1 = eta_ia[co.tor1];
		}
		
		if( co.tor2 == 0 ){
			tau2 = zero;
			eta2 = zero;
		} else {
			tau2 = tau_ia[co.tor2];
			eta2 = eta_ia[co.tor2];
		}

		Real const h = co.k/co.d2;

		// Iter over k,l torsions
		for( Size k = 1; k <= ntor(); ++k ){
			Vector eta = eta2[k] - eta1[k];
			Vector tau = tau2[k] - tau1[k];

			dk[k] = (co.c21).dot( eta ) - (co.r21).dot( tau );
			Real Hk = h*dk[k];

			/*
			printf("c21, r21: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
						 co.c21[0], co.c21[1], co.c21[2], co.r21[0], co.r21[1], co.r21[2] );

			printf("n,k %d %d, eta, tau, dk: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
						 i_co, k, eta[0], eta[1], eta[2],
						 tau[0], tau[1], tau[2],
						 dk[k] );
			*/

			for( Size l = 1; l <= k; ++l ){
				H[l][k] += Hk*dk[l];
			}
		}
	}

	// Fill upper diagonal
	for( Size k = 1; k <= ntor(); ++k ){
		for( Size l = 1; l <= k; ++l ){
			H[k][l] = H[l][k];
		}
	}

	TR.Debug << "Hessian" << std::endl;
	for( Size i = 1; i <= ntor(); ++i ){
		for( Size j = 1; j <= ntor(); ++j ){
			TR.Debug << " " << std::setw(9) << H[i][j];
		}
		TR.Debug << std::endl;
	}

	return H;
}

utility::vector1< utility::vector1< Real > >
NormalMode::convert_to_torsion_crd( 
		 utility::vector1< utility::vector1< Real > > const &U ){

	// First get COM coordinate
	Real Msum( 0.0 );
	Vector Rcom( 0.0, 0.0, 0.0 );
	for( Size iatm = 1; iatm <= natm(); ++iatm ){
		Msum += 1;
		Rcom += xyz( iatm );
	}
	Rcom /= Msum;

	// Translation to Center of Mass
	utility::vector1< Vector > xyz_com = xyz_;
	for( Size iatm = 1; iatm <= natm(); ++iatm ) xyz_com[iatm] -= Rcom;
	for( Size itor = 1; itor <= ntor(); ++itor ) s_[itor] -= Rcom;

	// Get Total inertia tensor: Assume mass is simply the same "1"
	numeric::MathMatrix< Real > Isum;
	numeric::MathMatrix< Real > Icorr_sum( 3, 3, 0.0 );

	for( Size iatm = 1; iatm <= natm(); ++iatm ){
		for( Size i = 0; i < 3; ++i )
			for( Size j = 0; j < 3; ++j )
				Icorr_sum[i][j] += xyz_com[iatm][i]*xyz_com[iatm][j];
	}
	// COM is at origin - Sum of coordinate should be zero
	Vector Ra_sum( 0.0, 0.0, 0.0 );
	Isum = calc_inertia_tensor( Icorr_sum, Ra_sum, Msum );

	/*
	std::cout << "Inertia tensor:" << std::endl;
	printf("%8.3f %8.3f %8.3f\n", Isum[0][0], Isum[0][1], Isum[0][2] );
	printf("%8.3f %8.3f %8.3f\n", Isum[1][0], Isum[1][1], Isum[1][2] );
	printf("%8.3f %8.3f %8.3f\n", Isum[2][0], Isum[2][1], Isum[2][2] );
	*/

	numeric::MathMatrix< Real > i_Isum = Isum.inverse_square_matrix();


	// Get Eckart Correction
	Real Ma( 0.0 ); // Mass upstream to a
	Vector Ra( 0.0, 0.0, 0.0 );  // Coordinate sum upstream to a
	Vector ta( 0.0, 0.0, 0.0 );  // Translation vector correction upstream to a
	Vector Aa( 0.0, 0.0, 0.0 );  // Rotation vector correction upstream to a
	numeric::MathMatrix< Real > Ia; // Inertia tensor upstream to a
	numeric::MathMatrix< Real > Icorr( 3, 3, 0.0 );
	numeric::MathMatrix< Real > J( 3*natm(), ntor(), 0.0 );
	numeric::MathMatrix< Real > Ut( 3*natm(), ntor(), 0.0 );

	for( Size a = ntor(); a >= 1; --a ){
		update_inertia_tensor( a, Ma, Ra, Icorr, Ia, xyz_com );

		// get Aa and ta
		if( eckart_correction_ ){
			calculate_Jacobi_correction( a, Ma, Ra, Msum, Ia, i_Isum, Aa, ta );
			// For debugging
			/*
				Vector local_rot = Aa + e_[a];
				Vector local_shift = ta - e_[a].cross( s_[a] );
				printf("a %d, local shift/rot global shift/rot: ", a);
				printf(" %8.3f %8.3f", std::sqrt( local_rot.dot_product( local_rot )), 
				std::sqrt(local_shift.dot_product( local_shift )) );
				printf(" %8.3f %8.3f", std::sqrt( Aa.dot( Aa ) ), std::sqrt( ta.dot( ta ) ) );
				printf("\n");
			*/
		}

		/// Generate Jacobian
		Vector torv( 0.0, 0.0, 0.0 );
		for( Size i = 1; i <= natm(); ++i ){
			if ( i > a_to_i( a ) ) torv = e_[a].cross( xyz_com[i] - s_[a] );

			// Aa == ta == 0 if eckart_correction_ == false
			Vector Jv = torv + Aa.cross( xyz_com[i] ) + ta; 

			// Be careful for index! J is 0 to n-1!
			J[3*i-3][a-1] = Jv[0];
			J[3*i-2][a-1] = Jv[1];
			J[3*i-1][a-1] = Jv[2];
		}
	}

	/*
	std::cout << "Jacobian" << std::endl;
	for( Size i = 0; i < ntor(); ++i ){
		for( Size j = 0; j < 3*natm(); ++j ){
			printf( "%8.3f ", J[j][i] );
		}
		printf( "\n" );
	}
	*/

	// Apply Jacobian to transform coordinate system ( Cartesian -> Torsion )
	// Convert into matrix system
	numeric::MathMatrix< Real > Utmp( 3*natm(), 3*natm(), 0.0 );
	for( Size i = 0; i < 3*natm(); ++i ){
		for( Size j = 0; j < 3*natm(); ++j ){
			Utmp[i][j] = U[i+1][j+1];
		}
	}

	// Convert Coordinate system
	Ut = J.transpose()*( Utmp*J );

	// Convert into vector1 system again
	utility::vector1< utility::vector1 < Real > > Ut1;
	Ut1.resize( ntor() );

	//std::cout << "Hessian, transformed" << std::endl;
	for( Size i = 0; i < ntor(); ++i ){
		Ut1[i+1].resize( ntor() );
		for( Size j = 0; j < ntor(); ++j ){
			Ut1[i+1][j+1] = Ut[i][j];
			//printf( "%8.3f ", Ut[i][j] );
		}
		//printf( "\n" );
	}

	return Ut1;
}

void 
NormalMode::update_inertia_tensor( Size const a, 
																	 Real &Ma,
																	 Vector &Ra, 
																	 numeric::MathMatrix< Real > &Icorr, 
																	 numeric::MathMatrix< Real > &Ia, 
																	 utility::vector1< Vector > const xyz_com 
																	 )
{
	Size start( a_to_i(a) + 1);
	Size end;
	if( a == ntor() ){
		end = natm();
	} else {
		end = a_to_i(a+1);
	}
	if( start <= end ){
		for( Size iatm = start; iatm <= end; ++iatm ){
			Ra += xyz_com[iatm];
			
			for( Size i = 0; i < 3; ++i)
				for( Size j = 0; j < 3; ++j)
					Icorr[i][j] += xyz_com[iatm][i]*xyz_com[iatm][j];

			Ma += 1.0;
			Ia = calc_inertia_tensor( Icorr, Ra, Ma );
		}
	}

	/*
	std::cout << "Inertia tensor:" << std::endl;
	printf("%8.3f %8.3f %8.3f\n", Ia[0][0], Ia[0][1], Ia[0][2] );
	printf("%8.3f %8.3f %8.3f\n", Ia[1][0], Ia[1][1], Ia[1][2] );
	printf("%8.3f %8.3f %8.3f\n", Ia[2][0], Ia[2][1], Ia[2][2] );
	*/
}

void 
NormalMode::calculate_Jacobi_correction( Size const a,
																				 Real const &Ma,
																				 Vector const &Ra, 
																				 Real const &Msum,
																				 numeric::MathMatrix< Real > const &Ia, 
																				 numeric::MathMatrix< Real > const &i_Isum,
																				 Vector &Aa,
																				 Vector &ta )
{

	// Rotation correction
	Vector IAa;
	Vector Rca( Ra/Ma );
	IAa[0] = -Ia[0][0]*e_[a][0] - Ia[0][1]*e_[a][1] - Ia[0][2]*e_[a][2];
	IAa[1] = -Ia[1][0]*e_[a][0] - Ia[1][1]*e_[a][1] - Ia[1][2]*e_[a][2];
	IAa[2] = -Ia[2][0]*e_[a][0] - Ia[2][1]*e_[a][1] - Ia[2][2]*e_[a][2];

	Vector const e_cross_Rca_sa( e_[a].cross( Rca - s_[a] ) );
	Vector tmpV( - Ra.cross( e_cross_Rca_sa ) );

	IAa += tmpV;

	Aa[0] = i_Isum[0][0]*IAa[0] + i_Isum[0][1]*IAa[1] + i_Isum[0][2]*IAa[2];
	Aa[1] = i_Isum[1][0]*IAa[0] + i_Isum[1][1]*IAa[1] + i_Isum[1][2]*IAa[2];
	Aa[2] = i_Isum[2][0]*IAa[0] + i_Isum[2][1]*IAa[1] + i_Isum[2][2]*IAa[2];

	// Translation correction
	//ta = -Aa.cross( Rcom ) - (Ma/Msum)*e_cross_Ra_sa - note that Rcom = O;
	ta = -(Ma/Msum)*e_cross_Rca_sa;

}

// Inertia tensor
// Assume mass is all "1"
numeric::MathMatrix<Real>
NormalMode::calc_inertia_tensor( numeric::MathMatrix<Real> const Icorr,
																 Vector const &Ra,
																 Real const Msum )
{

	numeric::MathMatrix<Real> Ia( 3, 3, 0.0 );

	for( Size i = 0; i < 3; ++i )
		for( Size j = 0; j < 3; ++j )
			Ia[i][j] = Icorr[i][j] - Ra[i]*Ra[j]/Msum;

	Real r2 = Ia[0][0] + Ia[1][1] + Ia[2][2];
	
	Ia[0][0] = r2 - Ia[0][0];
	Ia[1][1] = r2 - Ia[1][1];
	Ia[2][2] = r2 - Ia[2][2];
	Ia[0][1] = -Ia[0][1];
	Ia[0][2] = -Ia[0][2];
	Ia[1][2] = -Ia[1][2];
	Ia[1][0] = Ia[0][1];
	Ia[2][0] = Ia[0][2];
	Ia[2][1] = Ia[1][2]; 

	return Ia;
}

// compute (a2 + b2)^1/2 without destructive underflow or overflow
Real 
NormalMode::pythag(Real a, Real b)
{
  Real absa, absb;
  absa = std::fabs(a);
  absb = std::fabs(b);
  if (absa > absb) {
    return absa*std::sqrt(1.0+(absb/absa)*(absb/absa));
  } else {
    return (absb == 0.0 ? 0.0 : absb*std::sqrt(1.0+(absa/absb)*(absa/absb)));
  }
}

/******************************************************************************/
void 
NormalMode::svdcmp( utility::vector1< utility::vector1< Real > > &a, 
	    Size const m, Size const n, 
	    utility::vector1< Real > &w,
	    utility::vector1< utility::vector1< Real > > &v )
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
  Size flag,i,its,j,jj,k,l,nm;
  Real anorm,c,f,g,h,s,scale,x,y,z;
  
  //rv1 = dvector(1,n);
  utility::vector1< Real > rv1( n, 0.0 );
  g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */

  for (i=1; i<=n; i++) {

    l = i+1;
    rv1[i] = scale*g;

    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i <= m) {
      for (k=i;k<=m;k++){
				scale += std::fabs(a[k][i]);
      }
      if (scale) {
				for (k=i; k<=m; k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a[i][i] = f-g;
				for (j=l; j<=n; j++) {
					for (s=0.0,k=i; k<=m; k++){ s += a[k][i]*a[k][j]; }
					f=s/h;
					for (k=i; k<=m; k++){ a[k][j] += f*a[k][i]; }
				}
				for (k=i; k<=m; k++){ a[k][i] *= scale; }
      }
    }
		
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += std::fabs(a[i][k]);
      if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -sign(std::sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }

    anorm = std::max(anorm,(std::fabs(w[i])+std::fabs(rv1[i])));
  }

  for (i=n; i>=1; i--) { /* Accumulation of right-hand transformations. */
    if (i < n) {
      if (g) {
				for (j=l;j<=n;j++) /* Real division to avoid possible underflow. */
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
	
  for (i=std::min(m,n);i>=1;i--) { /* Accumulation of left-hand transformations. */
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
	
  for (k=n;k>=1;k--) { /* Diagonalization of the bidiagonal form. */
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) { /* Test for splitting. */
				nm=l-1; /* Note that rv1[1] is always zero. */
				if ((Real)(std::fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((Real)(std::fabs(w[nm])+anorm) == anorm) break;
      }
			
      if (flag) {
				c=0.0; /* Cancellation of rv1[l], if l > 1. */
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((Real)(std::fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
      }
			
      z=w[k];
      if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
      }
			
      x=w[l]; /* Shift from bottom 2-by-2 minor. */
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
      c=s=1.0; /* Next QR transformation: */
      for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
	
}

void
NormalMode::eigsrt( 
				utility::vector1< utility::vector1 < Real > > &eigvec,
				utility::vector1< Real > &eigval ){

	Size n( eigval.size() );
	Size k;

	for (Size i=1; i < n; i++) {
		Real p = eigval[k=i];

		for(Size j=i+1; j<=n; j++ ){
			if (eigval[j] <= p) p=eigval[k=j];
		}

		if (k != i) {
			eigval[k] = eigval[i];
			eigval[i] = p;

			for( Size j=1; j<=n; j++ ){
				// substitute
				Real tmp( eigvec[j][i] );
				eigvec[j][i] = eigvec[j][k];
				eigvec[j][k] = tmp;
			}
		}
	}
}

} // normalmode
} // protocols
