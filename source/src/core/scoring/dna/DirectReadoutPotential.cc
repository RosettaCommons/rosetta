// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/DirectReadoutPotential.cc
/// @brief  1st pass implementation of Kono + Sarai's protein-DNA interaction potential
/// @details  Needs polishing, converting to mini standards in some respects, but still in trial stage.
/// @author Amy Ticoll


// Unit Headers
#include <core/scoring/dna/DirectReadoutPotential.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>

// Package headers

// Project headers

#include <basic/database/open.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace dna {

/// @details ctor, reads data file. Need to configure to allow alternate tables/atom_sets
DirectReadoutPotential::DirectReadoutPotential()
:wt( 1.0/20 ),
	RT( 0.582 )
{
	fill_bins(A_bins, 'A');
	fill_bins(C_bins, 'C');
	fill_bins(G_bins, 'G');
	fill_bins(T_bins, 'T');
	get_pairs();

	aas_at_grid = 0;
	for ( int x=0; x<9; x++ ) {
		for ( int y=0; y<9; y++ ) {
			for ( int z=0; z<4; z++ ) {
				int num_aas =
					A_bins[x][y][z].length() +
					C_bins[x][y][z].length() +
					G_bins[x][y][z].length() +
					T_bins[x][y][z].length(); // need one letter code in list for this
				aas_at_grid += num_aas;
			}
		}
	}

	Real d_eab;
	for ( int x_bin=0; x_bin<9; x_bin++ ) {
		for ( int y_bin=0; y_bin<9; y_bin++ ) {
			for ( int z_bin=0; z_bin<4; z_bin++ ) {
				for ( int p=0; p<20; p++ ) {
					for ( int d=0; d<4; d++ ) {
						int num_aas =
							A_bins[x_bin][y_bin][z_bin].length() +
							C_bins[x_bin][y_bin][z_bin].length() +
							G_bins[x_bin][y_bin][z_bin].length() +
							T_bins[x_bin][y_bin][z_bin].length(); // need one letter code in list for this

						Real fs = Real(num_aas)/Real(aas_at_grid);
						int mab = num_pairs[d][p]; // m_pairs uses numbers, starting at 0

						int num_of_aa = 0;
						string lib_aa_list;
						if ( d==0 ) {  lib_aa_list=G_bins[x_bin][y_bin][z_bin];  } // follows order in AA.hh
						else if ( d==1 ) {  lib_aa_list=A_bins[x_bin][y_bin][z_bin];  }
						else if ( d==2 ) {  lib_aa_list=C_bins[x_bin][y_bin][z_bin];  }
						else {  lib_aa_list=T_bins[x_bin][y_bin][z_bin];  }

						for ( Size i=0; i<lib_aa_list.length(); i++ ) {
							int const aa( chemical::aa_from_oneletter_code( lib_aa_list[i] ) );
							if ( aa == p ) { num_of_aa++; } //***need to convert to_comp to its aa_number-1
						}

						Real fscore;
						Real gabs = Real(num_of_aa)/Real(mab);
						Real fabs = fs/(1 + mab*wt) + (mab*wt*gabs)/(1+mab*wt);
						if ( fs!=0 )  {  fscore = fabs/fs;  }
						else  {  fscore = 0;  }
						if ( fscore!=0 )  {  d_eab = -RT*log(fscore);  }
						else  {  d_eab = 0;  }

						score[x_bin][y_bin][z_bin][p][d]=d_eab;
					}
				}
			}
		}
	}
}

Real
DirectReadoutPotential::rsd_rsd_energy(conformation::Residue const & rsd1,conformation::Residue const & rsd2) const
{
	debug_assert( rsd1.is_protein() && rsd2.is_DNA() );

	// define coordinate frame
	bool const AG( rsd2.aa() == chemical::na_ade || rsd2.aa() == chemical::na_gua );
	Vector const origin( AG ? rsd2.xyz("N9") : rsd2.xyz("N1") );
	Vector p( AG ? rsd2.xyz("C4") : rsd2.xyz("C2") );
	Vector const x( ( p - origin ).normalized() );
	Vector z,y;
	z = ( dna::get_z_axis( rsd2, x ) );
	y = ( z.cross( x ) );
	Vector const calpha( rsd1.xyz("CA") - origin );
	Real const xx = dot(calpha,x);
	Real const yy = dot(calpha,y);
	Real const zz = dot(calpha,z);

	if ( ( xx >= -13.5 && xx <= 13.5  ) &&
			( yy >= -13.5 && yy <= 13.5  ) &&
			( zz >=  -6.0 && zz <=  6.0  ) ) {
		int const aa_bin = rsd1.aa() - 1;
		int const na_bin = rsd2.aa() - chemical::first_DNA_aa;

		int x_bin = get_xy_bin(xx);
		int y_bin = get_xy_bin(yy);
		int z_bin = get_z_bin(zz);
		debug_assert( aa_bin >= 0 && aa_bin < 20 && na_bin >= 0 && na_bin < 4 );
		return score[x_bin][y_bin][z_bin][aa_bin][na_bin];
	} else {
		return 0.0;
	}

}


void
DirectReadoutPotential::fill_bins(string (&my_array)[9][9][4], char const base )
{
	utility::io::izstream myfile;

	// open the file
	std::string const bins_filename( std::string("scoring/dna/") + base + "_bins.txt" );
	basic::database::open( myfile, bins_filename );

	if ( !myfile.good() ) utility_exit_with_message( "Unable to open file: "+bins_filename );

	for ( int i=0; i<9; ++i ) {
		for ( int j=0; j<9; ++j ) {
			for ( int k=0; k<4; ++k ) {
				std::string line, aas;
				int x,y,z;
				getline( myfile, line );
				std::istringstream l( line );
				l >> x >> y >> z >> aas;
				debug_assert( x == i && y == j && z == k );
				if ( aas == "-" ) {
					my_array[x][y][z]= "";
				} else {
					my_array[x][y][z]= aas;
				}
			}
		}
	}
	myfile.close();
}


void
DirectReadoutPotential::get_pairs()
{
	utility::io::izstream myfile;
	basic::database::open( myfile, "scoring/dna/m_pairs.txt" );

	if ( !myfile.good() ) utility_exit_with_message( "Unable to open m_pairs.txt" );

	for ( Size i=0; i<4; ++i ) {
		for ( Size j=0; j<20; ++j ) {
			string line;
			getline( myfile, line );
			std::istringstream l(line );
			Size itmp, jtmp;
			l >> itmp >> jtmp >> num_pairs[i][j];
			debug_assert( itmp == i && jtmp == j && !l.fail() );
		}
	}
	myfile.close();
}

int
DirectReadoutPotential::get_xy_bin(Real coord) const
{
	if ( coord <= -13.5 || coord >= 13.5 )  {  return -99;  }
	Real c = coord - ( -13.5 );
	int bin = int(floor( c/3 ));
	return bin;
}

int
DirectReadoutPotential::get_z_bin(Real coord) const
{
	if ( coord <= -6 || coord >= 6 )  {  return -99;  }
	Real c = coord - ( -6 );
	int bin = int(floor( c/3 ));
	return bin;
}

// int
// DirectReadoutPotential::num_pairs(int b, int a)
// {
//  for (int i=0; i<80; i++)
//   {  if (pair_base[i]==b && pair_aa[i]==a)  {  return pair_num[i];  }
//  }
// }


} // ns dna
} // ns scoring
} // ns core
