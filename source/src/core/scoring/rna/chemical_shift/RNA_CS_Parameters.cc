// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file   core/scoring/rna/chemical_shift/RNA_CS_Parameters.cc
/// @brief
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#include <core/scoring/rna/chemical_shift/RNA_CS_Parameters.hh>
#include <ObjexxFCL/format.hh>

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace ObjexxFCL;

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {

/////////////////////////////////////////////////////////RNA_CS_residue_parameters Class////////////////////////////////////////////////////////////////////


RNA_CS_residue_parameters::RNA_CS_residue_parameters( chemical::AA const & res_aa ):
	res_aa_( res_aa ),
	maxatoms_( 40 ),
	// Numerical data which is not specific for one of the bases
	RCCO_( 2.13 ),  // RCCO Constant in the formula for calculating the ring current effect
	MACQ_( 5.368 ), // MACQ Constant for the q factor of the Magnetic Anisotropy (this is the paramagnetic component coefficient)
	MACR_( 1.967 )  // MACR Constant for the r factor of the Magnetic Anisotropy (this is the diamagnetic component coefficient)
	//Currently no including charge parameters:
	// CHAC a coefficent for charge calculation
	// CHBC b coefficent for charge calculation
	// CONV conversion factor in charge calculation
	// CHD0 dielectric constant 0th order
	// CHD1 dielectric constant 1st order
{

	// Standard configuration file for the cs calculation in DNA
	// Numbers for RC and MA parameters from F.Ribas Prado & C.Giessner Prettre
	// Parameters for the calculation of the ring current and atomic magnetic anisotryopy
	// contributions to magnetic shielding constants:Nucleic Acid Bases and Intercalating
	// agents:  THEOCHEM / Journal of Molecular Structure 76 (1981) pp 81-92
	// OSHI corrections: Jenny Cromsigt et al. J.Biomol.NMR (2001)

	utility::vector1< core::Real > indiv_real_atom_data( last_atomdesc, 0.0 ); //Important to set to 0.0!

	realatomdata_.clear();
	realatomdata_.assign( maxatoms_, indiv_real_atom_data );

	ring_intensity_.clear(); //RCI ring current intensity (relative to benzene?)
	ring_radius_.clear();    //RCR: Radius of the ring (Angstrom)
	ring_height_.clear();     //RCH: Distance of the ring current loops to the molecular plane (Angstrom)

	atomnames_.clear();

	//P, OP2 and OP1 are considered sugar atom by NUCHEMICS [rename to backbone?]

	if ( res_aa_ == chemical::na_rad ) {

		BASE_ = "ADE"; //NUCHEMICS abbreviation.

		num_rings_ = 2; //RCCA 2

		ring_intensity_.push_back( 0.9000 );  //RCI1 0.9000
		ring_intensity_.push_back( 0.6600 );  //RCI2 0.6600

		ring_radius_.push_back( 1.3430 ); //RCR1 1.3430
		ring_radius_.push_back( 1.1540 ); //RCR2 1.1540

		ring_height_.push_back( 0.5660 ); //RCH1 0.5660
		ring_height_.push_back( 0.5660 ); //RCH2 0.5660

		//      1     2      3      4     5      6     7     8     9      10
		//ATOM  N1    C2     N3     C4    C5     C6    N6    N7    C8     N9
		//XDIR  0     2      0      1     0      0     0     0     0      0
		//YDIR  0     2      1      2     0      0     0     0     0      0

		//RCL1  1     1      1      1     1      1     0     0     0      0
		//RCL2  0     0      0      1     1      0     0     1     1      1

		//MACA  1     1      1      1     1      1     1     1     1      1
		//MAQX  2.943 3.050  3.048  3.148 3.147  2.386 2.567 1.972 3.301  2.139
		//MAQW -0.524 0.404 -0.093 -0.069 0.042 -0.075 0.000 0.186 0.131 -0.051
		//MAQY  2.644 2.700  3.239  2.959 2.656  3.158 2.000 3.206 2.433  2.957
		//MAQZ  2.000 2.000  2.000  2.000 2.000  2.000 2.000 2.000 2.000  2.000
		//MARX  5.020 5.308  4.427  5.615 5.842  5.349 5.267 5.506 5.537  4.942
		//MARY  5.283 5.308  4.902  5.615 5.842  5.349 5.267 4.980 5.537  4.942
		//MARZ  5.206 5.494  5.212  5.648 5.716  5.515 4.606 5.252 5.609  4.443

		atomnames_.push_back( "N1" ); //1
		atomnames_.push_back( "C2" ); //2
		atomnames_.push_back( "N3" ); //3
		atomnames_.push_back( "C4" ); //4
		atomnames_.push_back( "C5" ); //5
		atomnames_.push_back( "C6" ); //6
		atomnames_.push_back( "N6" ); //7
		atomnames_.push_back( "N7" ); //8
		atomnames_.push_back( "C8" ); //9
		atomnames_.push_back( "N9" ); //10

		realatomdata_[1 ][rcl1] = 1.0;
		realatomdata_[2 ][xdir] = 2.0;    realatomdata_[2 ][ydir] = 2.0;    realatomdata_[2 ][rcl1] = 1.0;
		realatomdata_[3 ][ydir] = 1.0;    realatomdata_[3 ][rcl1] = 1.0;
		realatomdata_[4 ][xdir] = 1.0;    realatomdata_[4 ][ydir] = 2.0;    realatomdata_[4 ][rcl1] = 1.0;    realatomdata_[4 ][rcl2] = 1.0;
		realatomdata_[5 ][rcl1] = 1.0;    realatomdata_[5 ][rcl2] = 1.0;
		realatomdata_[6 ][rcl1] = 1.0;

		realatomdata_[8 ][rcl2] = 1.0;
		realatomdata_[9 ][rcl2] = 1.0;
		realatomdata_[10][rcl2] = 1.0;

		realatomdata_[1 ][maca] = 1.0;   realatomdata_[1 ][maqx] = 2.943; realatomdata_[1 ][maqw] =  - 0.524;  realatomdata_[1 ][maqy] = 2.644;
		realatomdata_[2 ][maca] = 1.0;   realatomdata_[2 ][maqx] = 3.050; realatomdata_[2 ][maqw] = 0.404;  realatomdata_[2 ][maqy] = 2.700;
		realatomdata_[3 ][maca] = 1.0;   realatomdata_[3 ][maqx] = 3.048; realatomdata_[3 ][maqw] =  - 0.093;  realatomdata_[3 ][maqy] = 3.239;
		realatomdata_[4 ][maca] = 1.0;   realatomdata_[4 ][maqx] = 3.148; realatomdata_[4 ][maqw] =  - 0.069;  realatomdata_[4 ][maqy] = 2.959;
		realatomdata_[5 ][maca] = 1.0;   realatomdata_[5 ][maqx] = 3.147; realatomdata_[5 ][maqw] = 0.042;  realatomdata_[5 ][maqy] = 2.656;
		realatomdata_[6 ][maca] = 1.0;   realatomdata_[6 ][maqx] = 2.386; realatomdata_[6 ][maqw] =  - 0.075;  realatomdata_[6 ][maqy] = 3.158;
		realatomdata_[7 ][maca] = 1.0;   realatomdata_[7 ][maqx] = 2.567; realatomdata_[7 ][maqw] = 0.000;  realatomdata_[7 ][maqy] = 2.000;
		realatomdata_[8 ][maca] = 1.0;   realatomdata_[8 ][maqx] = 1.972; realatomdata_[8 ][maqw] = 0.186;  realatomdata_[8 ][maqy] = 3.206;
		realatomdata_[9 ][maca] = 1.0;   realatomdata_[9 ][maqx] = 3.301; realatomdata_[9 ][maqw] = 0.131;  realatomdata_[9 ][maqy] = 2.433;
		realatomdata_[10][maca] = 1.0;   realatomdata_[10][maqx] = 2.139; realatomdata_[10][maqw] =  - 0.051;  realatomdata_[10][maqy] = 2.957;

		realatomdata_[1 ][maqz] = 2.000; realatomdata_[1 ][marx] = 5.020; realatomdata_[1 ][mary] = 5.283;   realatomdata_[1 ][marz] = 5.206;
		realatomdata_[2 ][maqz] = 2.000; realatomdata_[2 ][marx] = 5.308; realatomdata_[2 ][mary] = 5.308;   realatomdata_[2 ][marz] = 5.494;
		realatomdata_[3 ][maqz] = 2.000; realatomdata_[3 ][marx] = 4.427; realatomdata_[3 ][mary] = 4.902;   realatomdata_[3 ][marz] = 5.212;
		realatomdata_[4 ][maqz] = 2.000; realatomdata_[4 ][marx] = 5.615;  realatomdata_[4 ][mary] = 5.615;   realatomdata_[4 ][marz] = 5.648;
		realatomdata_[5 ][maqz] = 2.000; realatomdata_[5 ][marx] = 5.842; realatomdata_[5 ][mary] = 5.842;   realatomdata_[5 ][marz] = 5.716;
		realatomdata_[6 ][maqz] = 2.000; realatomdata_[6 ][marx] = 5.349;  realatomdata_[6 ][mary] = 5.349;   realatomdata_[6 ][marz] = 5.515;
		realatomdata_[7 ][maqz] = 2.000; realatomdata_[7 ][marx] = 5.267;  realatomdata_[7 ][mary] = 5.267;   realatomdata_[7 ][marz] = 4.606;
		realatomdata_[8 ][maqz] = 2.000; realatomdata_[8 ][marx] = 5.506;  realatomdata_[8 ][mary] = 4.980;   realatomdata_[8 ][marz] = 5.252;
		realatomdata_[9 ][maqz] = 2.000; realatomdata_[9 ][marx] = 5.537; realatomdata_[9 ][mary] = 5.537;   realatomdata_[9 ][marz] = 5.609;
		realatomdata_[10][maqz] = 2.000; realatomdata_[10][marx] = 4.942; realatomdata_[10][mary] = 4.942;   realatomdata_[10][marz] = 4.443;


		//     11     12     13    14     15    16     17    18    19     20    21  22
		//ATOM P      OP2    OP1   C1'    C2'   C3'    O3'   C4'   O4'    C5'   O5' O2'
		//SUGA 1      1      1     1      1     1      1     1     1      1     1   1
		atomnames_.push_back( "P"  );  realatomdata_[11][suga] = 1.0; //11
		atomnames_.push_back( "OP2" );  realatomdata_[12][suga] = 1.0; //12
		atomnames_.push_back( "OP1" );  realatomdata_[13][suga] = 1.0; //13
		atomnames_.push_back( "C1'" );  realatomdata_[14][suga] = 1.0; //14
		atomnames_.push_back( "C2'" );  realatomdata_[15][suga] = 1.0; //15
		atomnames_.push_back( "C3'" );  realatomdata_[16][suga] = 1.0; //16
		atomnames_.push_back( "O3'" );  realatomdata_[17][suga] = 1.0; //17
		atomnames_.push_back( "C4'" );  realatomdata_[18][suga] = 1.0; //18
		atomnames_.push_back( "O4'" );  realatomdata_[19][suga] = 1.0; //19
		atomnames_.push_back( "C5'" );  realatomdata_[20][suga] = 1.0; //20
		atomnames_.push_back( "O5'" );  realatomdata_[21][suga] = 1.0; //21
		atomnames_.push_back( "O2'" );  realatomdata_[22][suga] = 1.0; //22


		//     23    24    25    26    27    28    29    30    31   32   33
		//ATOM H1'   H2'   HO2'  H3'   H4'   H5'   H5''  H2    H61  H62  H8
		//CSCA 1     1     0     1     1     1     1     1     0    0    1     // CSCA has to be nonzero if the chemical shift must be calculated for the atom
		//SUGA 1     1     1     1     1     1     1     0     0    0    0     // SUGA has to be true if an atom counts to the sugar part
		//OSHI 5.02  4.66  0     4.62  4.35  4.38  4.09  8.43  0    0    8.64  // OSHI, offset of the chemical shift calculation

		atomnames_.push_back( "H1'" );  realatomdata_[23][csca] = 1.0;  realatomdata_[23][suga] = 1.0; realatomdata_[23][oshi] = 5.02;  //23
		atomnames_.push_back( "H2'" );  realatomdata_[24][csca] = 1.0;  realatomdata_[24][suga] = 1.0; realatomdata_[24][oshi] = 4.66;  //24
		atomnames_.push_back( "HO2'" );  realatomdata_[25][csca] = 0.0;  realatomdata_[25][suga] = 1.0;  realatomdata_[25][oshi] = 0.00;  //25
		atomnames_.push_back( "H3'" );  realatomdata_[26][csca] = 1.0;  realatomdata_[26][suga] = 1.0;  realatomdata_[26][oshi] = 4.62;  //26
		atomnames_.push_back( "H4'" );  realatomdata_[27][csca] = 1.0;  realatomdata_[27][suga] = 1.0;  realatomdata_[27][oshi] = 4.35;  //27
		atomnames_.push_back( "H5'" );  realatomdata_[28][csca] = 1.0;  realatomdata_[28][suga] = 1.0; realatomdata_[28][oshi] = 4.38;  //28
		atomnames_.push_back( "H5''" );  realatomdata_[29][csca] = 1.0;  realatomdata_[29][suga] = 1.0; realatomdata_[29][oshi] = 4.09;  //29
		atomnames_.push_back( "H2"  );  realatomdata_[30][csca] = 1.0;  realatomdata_[30][suga] = 0.0; realatomdata_[30][oshi] = 8.43;  //30
		atomnames_.push_back( "H61" );  realatomdata_[31][csca] = 0.0;  realatomdata_[31][suga] = 0.0; realatomdata_[31][oshi] = 0.00;  //31
		atomnames_.push_back( "H62" );  realatomdata_[32][csca] = 0.0;   realatomdata_[32][suga] = 0.0; realatomdata_[32][oshi] = 0.00;  //32
		atomnames_.push_back( "H8"  );  realatomdata_[33][csca] = 1.0;  realatomdata_[33][suga] = 0.0; realatomdata_[33][oshi] = 8.64;  //33


	} else if ( res_aa_ == chemical::na_rgu ) {

		BASE_ = "GUA"; //NUCHEMICS abbreviation.

		num_rings_ = 2; //RCCA 2

		ring_intensity_.push_back( 0.3000 );  //RCI1 0.3000
		ring_intensity_.push_back( 0.6350 );  //RCI2 0.6350

		ring_radius_.push_back( 1.3610 ); //RCR1 1.3610
		ring_radius_.push_back( 1.1540 ); //RCR2 1.1540

		ring_height_.push_back( 0.5660 ); //RCH1 0.5660
		ring_height_.push_back( 0.5660 ); //RCH2 0.5660

		//      1     2      3      4     5      6     7     8     9      10     11
		//ATOM  N1    C2     N2    N3    C4     C5     C6    O6    N7     C8     N9
		//XDIR  0     2      0     0     1      0      0     0     0      0      0
		//YDIR  0     2      0     1     2      0      0     0     0      0      0
		//RCL1  1     1      0     1     1      1      1     0     0      0      0
		//RCL2  0     0      0     0     1      1      0     0     1      1      1
		//MACA  1     1      1     1     1      1      1     1     1      1      1
		//MAQX  2.731 2.985  2.135 1.862 3.200  3.194  3.327 2.809 2.054  3.289  2.330
		//MAQW -0.209 0.182 -0.234 0.181 0.025 -0.096 -0.014 0.000 0.141 -0.122 -0.051
		//MAQY  2.369 3.126  2.405 3.129 2.859  2.843  2.700 1.144 3.255  2.429  2.990
		//MAQZ  2.000 2.000  2.000 2.000 2.000  2.000  2.000 1.333 2.000  2.000  2.000
		//MARX  5.060 5.415  5.289 5.606 5.571  5.635  5.192 4.360 5.427  5.619  4.914
		//MARY  5.060 5.415  5.289 5.080 5.571  5.635  5.192 4.747 5.902  5.619  4.914
		//MARZ  4.503 5.548  4.617 5.301 5.716  5.758  5.436 4.692 5.212  5.650  4.429

		atomnames_.push_back( "N1" ); //1
		atomnames_.push_back( "C2" ); //2
		atomnames_.push_back( "N2" ); //3
		atomnames_.push_back( "N3" ); //4
		atomnames_.push_back( "C4" ); //5
		atomnames_.push_back( "C5" ); //6
		atomnames_.push_back( "C6" ); //7
		atomnames_.push_back( "O6" ); //8
		atomnames_.push_back( "N7" ); //9
		atomnames_.push_back( "C8" ); //10
		atomnames_.push_back( "N9" ); //11

		realatomdata_[1 ][rcl1] = 1.0;
		realatomdata_[2 ][xdir] = 2.0;    realatomdata_[2 ][ydir] = 2.0;    realatomdata_[2 ][rcl1] = 1.0;

		realatomdata_[4 ][ydir] = 1.0;    realatomdata_[4 ][rcl1] = 1.0;
		realatomdata_[5 ][xdir] = 1.0;    realatomdata_[5 ][ydir] = 2.0;    realatomdata_[5 ][rcl1] = 1.0;    realatomdata_[5 ][rcl2] = 1.0;
		realatomdata_[6 ][rcl1] = 1.0;    realatomdata_[6 ][rcl2] = 1.0;
		realatomdata_[7 ][rcl1] = 1.0;

		realatomdata_[9 ][rcl2] = 1.0;
		realatomdata_[10][rcl2] = 1.0;
		realatomdata_[11][rcl2] = 1.0;

		realatomdata_[1 ][maca] = 1.0;   realatomdata_[1 ][maqx] = 2.731; realatomdata_[1 ][maqw] =  - 0.209;  realatomdata_[1 ][maqy] = 2.369;
		realatomdata_[2 ][maca] = 1.0;   realatomdata_[2 ][maqx] = 2.985; realatomdata_[2 ][maqw] = 0.182;  realatomdata_[2 ][maqy] = 3.126;
		realatomdata_[3 ][maca] = 1.0;   realatomdata_[3 ][maqx] = 2.135; realatomdata_[3 ][maqw] =  - 0.234;  realatomdata_[3 ][maqy] = 2.405;
		realatomdata_[4 ][maca] = 1.0;   realatomdata_[4 ][maqx] = 1.862; realatomdata_[4 ][maqw] = 0.181;  realatomdata_[4 ][maqy] = 3.129;
		realatomdata_[5 ][maca] = 1.0;   realatomdata_[5 ][maqx] = 3.200; realatomdata_[5 ][maqw] = 0.025;  realatomdata_[5 ][maqy] = 2.859;
		realatomdata_[6 ][maca] = 1.0;   realatomdata_[6 ][maqx] = 3.194; realatomdata_[6 ][maqw] =  - 0.096;  realatomdata_[6 ][maqy] = 2.843;
		realatomdata_[7 ][maca] = 1.0;   realatomdata_[7 ][maqx] = 3.327; realatomdata_[7 ][maqw] =  - 0.014;  realatomdata_[7 ][maqy] = 2.700;
		realatomdata_[8 ][maca] = 1.0;   realatomdata_[8 ][maqx] = 2.809; realatomdata_[8 ][maqw] = 0.000;  realatomdata_[8 ][maqy] = 1.144;
		realatomdata_[9 ][maca] = 1.0;   realatomdata_[9 ][maqx] = 2.054; realatomdata_[9 ][maqw] = 0.141;  realatomdata_[9 ][maqy] = 3.255;
		realatomdata_[10][maca] = 1.0;   realatomdata_[10][maqx] = 3.289; realatomdata_[10][maqw] =  - 0.122;  realatomdata_[10][maqy] = 2.429;
		realatomdata_[11][maca] = 1.0;   realatomdata_[11][maqx] = 2.330; realatomdata_[11][maqw] =  - 0.051;  realatomdata_[11][maqy] = 2.990;


		realatomdata_[1 ][maqz] = 2.000; realatomdata_[1 ][marx] = 5.060; realatomdata_[1 ][mary] = 5.060;   realatomdata_[1 ][marz] = 4.503;
		realatomdata_[2 ][maqz] = 2.000; realatomdata_[2 ][marx] = 5.415; realatomdata_[2 ][mary] = 5.415;   realatomdata_[2 ][marz] = 5.548;
		realatomdata_[3 ][maqz] = 2.000; realatomdata_[3 ][marx] = 5.289; realatomdata_[3 ][mary] = 5.289;   realatomdata_[3 ][marz] = 4.617;
		realatomdata_[4 ][maqz] = 2.000; realatomdata_[4 ][marx] = 5.606;  realatomdata_[4 ][mary] = 5.080;   realatomdata_[4 ][marz] = 5.301;
		realatomdata_[5 ][maqz] = 2.000; realatomdata_[5 ][marx] = 5.571; realatomdata_[5 ][mary] = 5.571;   realatomdata_[5 ][marz] = 5.716;
		realatomdata_[6 ][maqz] = 2.000; realatomdata_[6 ][marx] = 5.635;  realatomdata_[6 ][mary] = 5.635;   realatomdata_[6 ][marz] = 5.758;
		realatomdata_[7 ][maqz] = 2.000; realatomdata_[7 ][marx] = 5.192;  realatomdata_[7 ][mary] = 5.192;   realatomdata_[7 ][marz] = 5.436;
		realatomdata_[8 ][maqz] = 1.333; realatomdata_[8 ][marx] = 4.360;  realatomdata_[8 ][mary] = 4.747;   realatomdata_[8 ][marz] = 4.692;
		realatomdata_[9 ][maqz] = 2.000; realatomdata_[9 ][marx] = 5.427; realatomdata_[9 ][mary] = 5.902;   realatomdata_[9 ][marz] = 5.212;
		realatomdata_[10][maqz] = 2.000; realatomdata_[10][marx] = 5.619; realatomdata_[10][mary] = 5.619;   realatomdata_[10][marz] = 5.650;
		realatomdata_[11][maqz] = 2.000; realatomdata_[11][marx] = 4.914; realatomdata_[11][mary] = 4.914;   realatomdata_[11][marz] = 4.429;

		//     12     13     14    15     16    17     18    19    20     21    22  23
		//ATOM P      OP2    OP1   C1'    C2'   C3'    O3'   C4'   O4'    C5'   O5' O2'
		//SUGA 1      1      1     1      1     1      1     1     1      1     1   1
		atomnames_.push_back( "P"  );  realatomdata_[12][suga] = 1.0; //12
		atomnames_.push_back( "OP2" );  realatomdata_[13][suga] = 1.0; //13
		atomnames_.push_back( "OP1" );  realatomdata_[14][suga] = 1.0; //14
		atomnames_.push_back( "C1'" );  realatomdata_[15][suga] = 1.0; //15
		atomnames_.push_back( "C2'" );  realatomdata_[16][suga] = 1.0; //16
		atomnames_.push_back( "C3'" );  realatomdata_[17][suga] = 1.0; //17
		atomnames_.push_back( "O3'" );  realatomdata_[18][suga] = 1.0; //18
		atomnames_.push_back( "C4'" );  realatomdata_[19][suga] = 1.0; //19
		atomnames_.push_back( "O4'" );  realatomdata_[20][suga] = 1.0; //20
		atomnames_.push_back( "C5'" );  realatomdata_[21][suga] = 1.0; //21
		atomnames_.push_back( "O5'" );  realatomdata_[22][suga] = 1.0; //22
		atomnames_.push_back( "O2'" );  realatomdata_[23][suga] = 1.0; //23


		//     24    25    26    27    28    29    30    31    32   33  34
		//ATOM H1'   H2'   HO2'  H3'   H4'   H5'   H5''  H1   H21  H22  H8
		//CSCA 1     1     0     1     1     1     1     0    0    0    1
		//SUGA 1     1     1     1     1     1     1     0    0    0    0
		//OSHI 5.20  4.66  0     4.62  4.35  4.38  4.09  0    0    0    8.10

		atomnames_.push_back( "H1'" );  realatomdata_[24][csca] = 1.0;  realatomdata_[24][suga] = 1.0; realatomdata_[24][oshi] = 5.20;  //24
		atomnames_.push_back( "H2'" );  realatomdata_[25][csca] = 1.0;  realatomdata_[25][suga] = 1.0; realatomdata_[25][oshi] = 4.66;  //25
		atomnames_.push_back( "HO2'" );  realatomdata_[26][csca] = 0.0;  realatomdata_[26][suga] = 1.0;  realatomdata_[26][oshi] = 0.00;  //26
		atomnames_.push_back( "H3'" );  realatomdata_[27][csca] = 1.0;  realatomdata_[27][suga] = 1.0;  realatomdata_[27][oshi] = 4.62;  //27
		atomnames_.push_back( "H4'" );  realatomdata_[28][csca] = 1.0;  realatomdata_[28][suga] = 1.0;  realatomdata_[28][oshi] = 4.35;  //28
		atomnames_.push_back( "H5'" );  realatomdata_[29][csca] = 1.0;  realatomdata_[29][suga] = 1.0; realatomdata_[29][oshi] = 4.38;  //29
		atomnames_.push_back( "H5''" );  realatomdata_[30][csca] = 1.0;  realatomdata_[30][suga] = 1.0; realatomdata_[30][oshi] = 4.09;  //30
		atomnames_.push_back( "H1"  );  realatomdata_[31][csca] = 0.0;  realatomdata_[31][suga] = 0.0; realatomdata_[31][oshi] = 0.00;  //31
		atomnames_.push_back( "H22" );  realatomdata_[32][csca] = 0.0;  realatomdata_[32][suga] = 0.0; realatomdata_[32][oshi] = 0.00;  //32
		atomnames_.push_back( "H21" );  realatomdata_[33][csca] = 0.0;  realatomdata_[33][suga] = 0.0; realatomdata_[33][oshi] = 0.00;  //33
		atomnames_.push_back( "H8"  );  realatomdata_[34][csca] = 1.0;  realatomdata_[34][suga] = 0.0; realatomdata_[34][oshi] = 8.10;  //34


	} else if ( res_aa_ == chemical::na_rcy ) {

		BASE_ = "CYT"; //BASE CYT //NUCHEMICS abbreviation.

		num_rings_ = 1; //RCCA 1

		ring_intensity_.push_back( 0.2750 );  //RCI1 0.2750

		ring_radius_.push_back( 1.3675 ); //RCR1 1.3675

		ring_height_.push_back( 0.5770 ); //RCH1 0.5770

		//      1      2      3      4      5     6     7      8
		//ATOM  N1     C2     O2     N3     C4    N4    C5     C6
		//XDIR  0      2      0      0      0     0     0      1
		//YDIR  1      2      0      0      0     0     0      2
		//RCL1  1      1      0      1      1     0     1      1
		//MACA  1      1      1      1      1     1     1      1
		//MAQX  2.675  2.917  1.496  2.275  3.284 2.625 3.240  2.998
		//MAQW -0.592 -0.236 -0.724 -0.032 -0.255 0.000 0.213 -0.221
		//MAQY  2.715  3.179  2.264  2.826  2.441 2.000 2.396  3.120
		//MAQZ  2.000  2.000  1.333  2.000  2.000 2.000 2.000  2.000
		//MARX  5.034  5.211  4.688  5.017  5.390 5.209 5.985  5.571
		//MARY  5.034  5.211  4.495  5.334  5.390 5.209 5.985  5.571
		//MARZ  4.489  5.445  4.711  5.231  5.535 4.577 5.833  5.626

		atomnames_.push_back( "N1" ); //1
		atomnames_.push_back( "C2" ); //2
		atomnames_.push_back( "O2" ); //3
		atomnames_.push_back( "N3" ); //4
		atomnames_.push_back( "C4" ); //5
		atomnames_.push_back( "N4" ); //6
		atomnames_.push_back( "C5" ); //7
		atomnames_.push_back( "C6" ); //8


		realatomdata_[1 ][ydir] = 1.0;    realatomdata_[1 ][rcl1] = 1.0;
		realatomdata_[2 ][xdir] = 2.0;    realatomdata_[2 ][ydir] = 2.0;    realatomdata_[2 ][rcl1] = 1.0;

		realatomdata_[4 ][rcl1] = 1.0;
		realatomdata_[5 ][rcl1] = 1.0;

		realatomdata_[7 ][rcl1] = 1.0;
		realatomdata_[8 ][xdir] = 1.0;    realatomdata_[8 ][ydir] = 2.0;    realatomdata_[8 ][rcl1] = 1.0;


		realatomdata_[1 ][maca] = 1.0;   realatomdata_[1 ][maqx] = 2.675; realatomdata_[1 ][maqw] =  - 0.592;  realatomdata_[1 ][maqy] = 2.715;
		realatomdata_[2 ][maca] = 1.0;   realatomdata_[2 ][maqx] = 2.917; realatomdata_[2 ][maqw] =  - 0.236;  realatomdata_[2 ][maqy] = 3.179;
		realatomdata_[3 ][maca] = 1.0;   realatomdata_[3 ][maqx] = 1.496; realatomdata_[3 ][maqw] =  - 0.724;  realatomdata_[3 ][maqy] = 2.264;
		realatomdata_[4 ][maca] = 1.0;   realatomdata_[4 ][maqx] = 2.275; realatomdata_[4 ][maqw] =  - 0.032;  realatomdata_[4 ][maqy] = 2.826;
		realatomdata_[5 ][maca] = 1.0;   realatomdata_[5 ][maqx] = 3.284; realatomdata_[5 ][maqw] =  - 0.255;  realatomdata_[5 ][maqy] = 2.441;
		realatomdata_[6 ][maca] = 1.0;   realatomdata_[6 ][maqx] = 2.625; realatomdata_[6 ][maqw] = 0.000;  realatomdata_[6 ][maqy] = 2.000;
		realatomdata_[7 ][maca] = 1.0;   realatomdata_[7 ][maqx] = 3.240; realatomdata_[7 ][maqw] = 0.213;  realatomdata_[7 ][maqy] = 2.396;
		realatomdata_[8 ][maca] = 1.0;   realatomdata_[8 ][maqx] = 2.998; realatomdata_[8 ][maqw] =  - 0.221;  realatomdata_[8 ][maqy] = 3.120;

		realatomdata_[1 ][maqz] = 2.000; realatomdata_[1 ][marx] = 5.034; realatomdata_[1 ][mary] = 5.034;   realatomdata_[1 ][marz] = 4.489;
		realatomdata_[2 ][maqz] = 2.000; realatomdata_[2 ][marx] = 5.211; realatomdata_[2 ][mary] = 5.211;   realatomdata_[2 ][marz] = 5.445;
		realatomdata_[3 ][maqz] = 1.333; realatomdata_[3 ][marx] = 4.688; realatomdata_[3 ][mary] = 4.495;   realatomdata_[3 ][marz] = 4.711;
		realatomdata_[4 ][maqz] = 2.000; realatomdata_[4 ][marx] = 5.017;  realatomdata_[4 ][mary] = 5.334;   realatomdata_[4 ][marz] = 5.231;
		realatomdata_[5 ][maqz] = 2.000; realatomdata_[5 ][marx] = 5.390; realatomdata_[5 ][mary] = 5.390;   realatomdata_[5 ][marz] = 5.535;
		realatomdata_[6 ][maqz] = 2.000; realatomdata_[6 ][marx] = 5.209;  realatomdata_[6 ][mary] = 5.209;   realatomdata_[6 ][marz] = 4.577;
		realatomdata_[7 ][maqz] = 2.000; realatomdata_[7 ][marx] = 5.985;  realatomdata_[7 ][mary] = 5.985;   realatomdata_[7 ][marz] = 5.833;
		realatomdata_[8 ][maqz] = 2.000; realatomdata_[8 ][marx] = 5.571;  realatomdata_[8 ][mary] = 5.571;   realatomdata_[8 ][marz] = 5.626;

		//     9      10     11    12     13    14     15    16    17    18    19   20
		//ATOM P      OP2    OP1   C1'    C2'   C3'    O3'   C4'   O4'   C5'   O5'  O2'
		//SUGA 1      1      1     1      1     1      1     1     1     1     1    1
		atomnames_.push_back( "P"  );  realatomdata_[9 ][suga] = 1.0; //9
		atomnames_.push_back( "OP2" );  realatomdata_[10][suga] = 1.0; //10
		atomnames_.push_back( "OP1" );  realatomdata_[11][suga] = 1.0; //11
		atomnames_.push_back( "C1'" );  realatomdata_[12][suga] = 1.0; //12
		atomnames_.push_back( "C2'" );  realatomdata_[13][suga] = 1.0; //13
		atomnames_.push_back( "C3'" );  realatomdata_[14][suga] = 1.0; //14
		atomnames_.push_back( "O3'" );  realatomdata_[15][suga] = 1.0; //15
		atomnames_.push_back( "C4'" );  realatomdata_[16][suga] = 1.0; //16
		atomnames_.push_back( "O4'" );  realatomdata_[17][suga] = 1.0; //17
		atomnames_.push_back( "C5'" );  realatomdata_[18][suga] = 1.0; //18
		atomnames_.push_back( "O5'" );  realatomdata_[19][suga] = 1.0; //19
		atomnames_.push_back( "O2'" );  realatomdata_[20][suga] = 1.0; //20

		//     21    22    23    24    25    26    27    28   29   30    31
		//ATOM H1'   H2'   HO2'  H3'   H4'   H5'   H5''  H41  H42  H5    H6
		//CSCA 1     1     0     1     1     1     1     0    0    1     1
		//SUGA 1     1     1     1     1     1     1     0    0    0     0
		//OSHI 5.23  4.66  0     4.62  4.35  4.38  4.09  0    0    6.12  8.16
		atomnames_.push_back( "H1'" );  realatomdata_[21][csca] = 1.0;  realatomdata_[21][suga] = 1.0; realatomdata_[21][oshi] = 5.23;  //21
		atomnames_.push_back( "H2'" );  realatomdata_[22][csca] = 1.0;  realatomdata_[22][suga] = 1.0; realatomdata_[22][oshi] = 4.66;  //22
		atomnames_.push_back( "HO2'" );  realatomdata_[23][csca] = 0.0;  realatomdata_[23][suga] = 1.0;  realatomdata_[23][oshi] = 0.00;  //23
		atomnames_.push_back( "H3'" );  realatomdata_[24][csca] = 1.0;  realatomdata_[24][suga] = 1.0;  realatomdata_[24][oshi] = 4.62;  //24
		atomnames_.push_back( "H4'" );  realatomdata_[25][csca] = 1.0;  realatomdata_[25][suga] = 1.0;  realatomdata_[25][oshi] = 4.35;  //25
		atomnames_.push_back( "H5'" );  realatomdata_[26][csca] = 1.0;  realatomdata_[26][suga] = 1.0; realatomdata_[26][oshi] = 4.38;  //26
		atomnames_.push_back( "H5''" );  realatomdata_[27][csca] = 1.0;  realatomdata_[27][suga] = 1.0; realatomdata_[27][oshi] = 4.09;  //27
		atomnames_.push_back( "H41" );  realatomdata_[28][csca] = 0.0;  realatomdata_[28][suga] = 0.0; realatomdata_[28][oshi] = 0.00;  //28
		atomnames_.push_back( "H42" );  realatomdata_[29][csca] = 0.0;  realatomdata_[29][suga] = 0.0; realatomdata_[29][oshi] = 0.00;  //29
		atomnames_.push_back( "H5"  );  realatomdata_[30][csca] = 1.0;  realatomdata_[30][suga] = 0.0; realatomdata_[30][oshi] = 6.12;  //30
		atomnames_.push_back( "H6"  );  realatomdata_[31][csca] = 1.0;  realatomdata_[31][suga] = 0.0; realatomdata_[31][oshi] = 8.16;  //31

	} else if ( res_aa_ == chemical::na_ura ) {


		BASE_ = "URI"; //BASE URI //NUCHEMICS abbreviation.

		num_rings_ = 1; //RCCA 1

		ring_intensity_.push_back( 0.1110 );  //RCI1 0.1110

		ring_radius_.push_back( 1.3790 ); //RCR1 1.3790

		ring_height_.push_back( 0.5770 ); //RCH1 0.5770

		//      1      2      3      4      5     6     7      8
		//ATOM  N1     C2     O2     N3     C4    O4    C5     C6
		//XDIR  0      2      0      0      0     0     0      1
		//YDIR  1      2      0      0      0     0     0      2
		//RCL1  1      1      0      1      1     0     1      1
		//MACA  1      1      1      1      1     1     1      1
		//MAQX  2.257  2.921  1.517  2.653  3.348 2.857 3.275  3.299
		//MAQW -0.020 -0.253 -0.726 -0.215 -0.024 0.000 0.192 -0.233
		//MAQY  2.773  3.176  2.269  2.373  2.705 1.232 2.332  2.404
		//MAQZ  2.000  2.000  1.333  2.000  2.000 1.333 2.000  2.000
		//MARX  5.097  5.240  4.677  5.112  5.260 4.309 5.883  5.646
		//MARY  5.097  5.240  4.484  5.112  5.260 4.696 5.883  5.646
		//MARZ  4.521  5.460  4.705  4.529  5.470 4.666 5.782  5.663

		atomnames_.push_back( "N1" ); //1
		atomnames_.push_back( "C2" ); //2
		atomnames_.push_back( "O2" ); //3
		atomnames_.push_back( "N3" ); //4
		atomnames_.push_back( "C4" ); //5
		atomnames_.push_back( "O4" ); //6
		atomnames_.push_back( "C5" ); //7
		atomnames_.push_back( "C6" ); //8


		realatomdata_[1 ][ydir] = 1.0;    realatomdata_[1 ][rcl1] = 1.0;
		realatomdata_[2 ][xdir] = 2.0;    realatomdata_[2 ][ydir] = 2.0;    realatomdata_[2 ][rcl1] = 1.0;

		realatomdata_[4 ][rcl1] = 1.0;
		realatomdata_[5 ][rcl1] = 1.0;

		realatomdata_[7 ][rcl1] = 1.0;
		realatomdata_[8 ][xdir] = 1.0;    realatomdata_[8 ][ydir] = 2.0;    realatomdata_[8 ][rcl1] = 1.0;


		realatomdata_[1 ][maca] = 1.0;   realatomdata_[1 ][maqx] = 2.257; realatomdata_[1 ][maqw] =  - 0.020;  realatomdata_[1 ][maqy] = 2.773;
		realatomdata_[2 ][maca] = 1.0;   realatomdata_[2 ][maqx] = 2.921; realatomdata_[2 ][maqw] =  - 0.253;  realatomdata_[2 ][maqy] = 3.176;
		realatomdata_[3 ][maca] = 1.0;   realatomdata_[3 ][maqx] = 1.517; realatomdata_[3 ][maqw] =  - 0.726;  realatomdata_[3 ][maqy] = 2.269;
		realatomdata_[4 ][maca] = 1.0;   realatomdata_[4 ][maqx] = 2.653; realatomdata_[4 ][maqw] =  - 0.215;  realatomdata_[4 ][maqy] = 2.373;
		realatomdata_[5 ][maca] = 1.0;   realatomdata_[5 ][maqx] = 3.348; realatomdata_[5 ][maqw] =  - 0.024;  realatomdata_[5 ][maqy] = 2.705;
		realatomdata_[6 ][maca] = 1.0;   realatomdata_[6 ][maqx] = 2.857; realatomdata_[6 ][maqw] = 0.000;  realatomdata_[6 ][maqy] = 1.232;
		realatomdata_[7 ][maca] = 1.0;   realatomdata_[7 ][maqx] = 3.275; realatomdata_[7 ][maqw] = 0.192;  realatomdata_[7 ][maqy] = 2.332;
		realatomdata_[8 ][maca] = 1.0;   realatomdata_[8 ][maqx] = 3.299; realatomdata_[8 ][maqw] =  - 0.233;  realatomdata_[8 ][maqy] = 2.404;

		realatomdata_[1 ][maqz] = 2.000; realatomdata_[1 ][marx] = 5.097; realatomdata_[1 ][mary] = 5.097;   realatomdata_[1 ][marz] = 4.521;
		realatomdata_[2 ][maqz] = 2.000; realatomdata_[2 ][marx] = 5.240; realatomdata_[2 ][mary] = 5.240;   realatomdata_[2 ][marz] = 5.460;
		realatomdata_[3 ][maqz] = 1.333; realatomdata_[3 ][marx] = 4.677; realatomdata_[3 ][mary] = 4.484;   realatomdata_[3 ][marz] = 4.705;
		realatomdata_[4 ][maqz] = 2.000; realatomdata_[4 ][marx] = 5.112;  realatomdata_[4 ][mary] = 5.112;   realatomdata_[4 ][marz] = 4.529;
		realatomdata_[5 ][maqz] = 2.000; realatomdata_[5 ][marx] = 5.260; realatomdata_[5 ][mary] = 5.260;   realatomdata_[5 ][marz] = 5.470;
		realatomdata_[6 ][maqz] = 1.333; realatomdata_[6 ][marx] = 4.309;  realatomdata_[6 ][mary] = 4.696;   realatomdata_[6 ][marz] = 4.666;
		realatomdata_[7 ][maqz] = 2.000; realatomdata_[7 ][marx] = 5.883;  realatomdata_[7 ][mary] = 5.883;   realatomdata_[7 ][marz] = 5.782;
		realatomdata_[8 ][maqz] = 2.000; realatomdata_[8 ][marx] = 5.646;  realatomdata_[8 ][mary] = 5.646;   realatomdata_[8 ][marz] = 5.663;

		//     9      10     11    12     13    14     15    16    17    18    19   20
		//ATOM P      OP2    OP1   C1'    C2'   C3'    O3'   C4'   O4'    C5'   O5' O2'
		//SUGA 1      1      1     1      1     1      1     1     1      1     1   1
		atomnames_.push_back( "P"  );  realatomdata_[9 ][suga] = 1.0; //9
		atomnames_.push_back( "OP2" );  realatomdata_[10][suga] = 1.0; //10
		atomnames_.push_back( "OP1" );  realatomdata_[11][suga] = 1.0; //11
		atomnames_.push_back( "C1'" );  realatomdata_[12][suga] = 1.0; //12
		atomnames_.push_back( "C2'" );  realatomdata_[13][suga] = 1.0; //13
		atomnames_.push_back( "C3'" );  realatomdata_[14][suga] = 1.0; //14
		atomnames_.push_back( "O3'" );  realatomdata_[15][suga] = 1.0; //15
		atomnames_.push_back( "C4'" );  realatomdata_[16][suga] = 1.0; //16
		atomnames_.push_back( "O4'" );  realatomdata_[17][suga] = 1.0; //17
		atomnames_.push_back( "C5'" );  realatomdata_[18][suga] = 1.0; //18
		atomnames_.push_back( "O5'" );  realatomdata_[19][suga] = 1.0; //19
		atomnames_.push_back( "O2'" );  realatomdata_[20][suga] = 1.0; //20

		//     21    22    23    24    25    26    27    28   29    30
		//ATOM H1'   H2'   HO2'  H3'   H4'   H5'   H5''  H3   H6    H5
		//CSCA 1     1     0     1     1     1     1     0    1     1
		//SUGA 1     1     1     1     1     1     1     0    0     0
		//OSHI 5.57  4.66  0     4.62  4.35  4.38  4.09  0    8.25  6.04
		atomnames_.push_back( "H1'" );  realatomdata_[21][csca] = 1.0;  realatomdata_[21][suga] = 1.0; realatomdata_[21][oshi] = 5.57;  //21
		atomnames_.push_back( "H2'" );  realatomdata_[22][csca] = 1.0;  realatomdata_[22][suga] = 1.0; realatomdata_[22][oshi] = 4.66;  //22
		atomnames_.push_back( "HO2'" );  realatomdata_[23][csca] = 0.0;  realatomdata_[23][suga] = 1.0;  realatomdata_[23][oshi] = 0.00;  //23
		atomnames_.push_back( "H3'" );  realatomdata_[24][csca] = 1.0;  realatomdata_[24][suga] = 1.0;  realatomdata_[24][oshi] = 4.62;  //24
		atomnames_.push_back( "H4'" );  realatomdata_[25][csca] = 1.0;  realatomdata_[25][suga] = 1.0;  realatomdata_[25][oshi] = 4.35;  //25
		atomnames_.push_back( "H5'" );  realatomdata_[26][csca] = 1.0;  realatomdata_[26][suga] = 1.0; realatomdata_[26][oshi] = 4.38;  //26
		atomnames_.push_back( "H5''" );  realatomdata_[27][csca] = 1.0;  realatomdata_[27][suga] = 1.0; realatomdata_[27][oshi] = 4.09;  //27
		atomnames_.push_back( "H3"  );  realatomdata_[28][csca] = 0.0;  realatomdata_[28][suga] = 0.0; realatomdata_[28][oshi] = 0.00;  //28
		atomnames_.push_back( "H6"  );  realatomdata_[29][csca] = 1.0;  realatomdata_[29][suga] = 0.0; realatomdata_[29][oshi] = 8.25;  //29
		atomnames_.push_back( "H5"  );  realatomdata_[30][csca] = 1.0;  realatomdata_[30][suga] = 0.0; realatomdata_[30][oshi] = 6.04;  //30


	} else {
		utility_exit_with_message( "Invalid res_aa_ ( " + ObjexxFCL::string_of( res_aa_ ) + " )!" );
	}

	///Consistency check////
	//if(num_rings_!=atoms_in_ring_list_.size()){
	// utility_exit_with_message("atom ("+ string_of(num_rings_) + ")>atomnames_.size() ("+ string_of(atoms_in_ring_list_.size()) + ")!");
	//}

	if ( num_rings_ != ring_intensity_.size() ) {
		utility_exit_with_message( "atom ( " + string_of( num_rings_ ) + " ) > ring_intensity_.size() ( " + string_of( ring_intensity_.size() ) + " )!" );
	}

	if ( num_rings_ != ring_radius_.size() ) {
		utility_exit_with_message( "atom ( " + string_of( num_rings_ ) + " ) > ring_radius_.size() ( " + string_of( ring_radius_.size() ) + " )!" );
	}

	if ( num_rings_ != ring_height_.size() ) {
		utility_exit_with_message( "atom ( " + string_of( num_rings_ ) + " ) > ring_height_.size() ( " + string_of( ring_height_.size() ) + " )!" );
	}

	if ( atomnames_.size() > maxatoms_ ) utility_exit_with_message( "atomnames_.size() > maxatoms_" );

	std::cout << " for BASE_ = " << BASE_ << " :";
	std::cout << " atomnames_.size() = " << atomnames_.size();
	std::cout << std::endl;
}


//destructor
RNA_CS_residue_parameters::~RNA_CS_residue_parameters(){}

////////////////////////////////////////////////////////////
std::string const
RNA_CS_residue_parameters::base_name() const{
	return BASE_;
}

////////////////////////////////////////////////////////////
Size
RNA_CS_residue_parameters::num_rings() const{
	return num_rings_;
}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::ring_intensity( Size const ring_ID ) const{

	if ( ring_ID > num_rings_ ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) > num_rings_( " + string_of( num_rings_ ) + " )!" );

	if ( ring_ID < 1 ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) < 1!" );


	return ring_intensity_[ring_ID];

}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::ring_radius( Size const ring_ID ) const{

	if ( ring_ID > num_rings_ ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) > num_rings_( " + string_of( num_rings_ ) + " )!" );

	if ( ring_ID < 1 ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) < 1!" );

	return ring_radius_[ring_ID];

}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::ring_height( Size const ring_ID ) const{

	if ( ring_ID > num_rings_ ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) > num_rings_( " + string_of( num_rings_ ) + " )!" );

	if ( ring_ID < 1 ) utility_exit_with_message( "ring_ID ( " + string_of( ring_ID ) + " ) < 1!" );

	return ring_height_[ring_ID];

}

////////////////////////////////////////////////////////////
Size
RNA_CS_residue_parameters::get_atomnames_size() const{

	return atomnames_.size();

}

////////////////////////////////////////////////////////////
std::string const
RNA_CS_residue_parameters::get_atomname( Size const count ) const{

	if ( count > atomnames_.size() ) utility_exit_with_message( "count ( " + string_of( count ) + " ) > atomnames_.size() ( " + string_of( atomnames_.size() ) + " )!" );

	return atomnames_[count];
}


////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::atom_data( Size const atom, atomitem const item ) const{


	if ( atom > atomnames_.size() ) utility_exit_with_message( "atom ( " + string_of( atom ) + " ) > atomnames_.size() ( " + string_of( atomnames_.size() ) + " )!" );

	if ( item > last_atomdesc ) utility_exit_with_message( "atom ( " + string_of( item ) + " ) > atomnames_.size() ( " + string_of( last_atomdesc ) + " )!" );

	if ( item < 1 ) utility_exit_with_message( "item ( " + ObjexxFCL::string_of( item ) + " ) < 1!" );


	return realatomdata_[atom][item];


}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::ring_current_coeff() const
{
	return RCCO_;

}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::magentic_anisotropy_r_coeff() const
{
	return MACR_;

}

////////////////////////////////////////////////////////////
Real
RNA_CS_residue_parameters::magentic_anisotropy_q_coeff() const
{
	return MACQ_;

}
////////////////////////////////////////////////////////////
chemical::AA
RNA_CS_residue_parameters::aa() const{
	return res_aa_;
}

/////////////////////////////////////////////////////////RNA_CS_parameters Class////////////////////////////////////////////////////////////////////

//constructor!
RNA_CS_parameters::RNA_CS_parameters():
	COMM_( "RNA data set without charge" ),
	//From datafile.cc..include these when needed.
	//IGNORE:const int maxdiffbases = 10; ##NUMBER OF TYPE OF BASES WHICH FOR RNA SHOULD BE FOUR (RAD, RGU, RCY, URA)
	//IGNORE:const int maxbasenamelength = 3;
	//IGNORE:const int maxatomnamelength = 5;
	//const int maxstrandnamelength = 6;
	//const int maxfilenamelength = 256;
	//const char remarkstart = '#';

	CS_RAD_params_( RNA_CS_residue_parameters( chemical::na_rad ) ),
	CS_RGU_params_( RNA_CS_residue_parameters( chemical::na_rgu ) ),
	CS_RCY_params_( RNA_CS_residue_parameters( chemical::na_rcy ) ),
	CS_URA_params_( RNA_CS_residue_parameters( chemical::na_ura ) )

{
}


//destructor
RNA_CS_parameters::~RNA_CS_parameters(){}

////////////////////////////////////////////////////////////
RNA_CS_residue_parameters const &
RNA_CS_parameters::get_RNA_CS_residue_parameters( chemical::AA const res_aa ) const
{

	if ( res_aa == chemical::na_rad ) return CS_RAD_params_;

	if ( res_aa == chemical::na_rgu ) return CS_RGU_params_;

	if ( res_aa == chemical::na_rcy ) return CS_RCY_params_;

	if ( res_aa == chemical::na_ura ) return CS_URA_params_;

	utility_exit_with_message( "Invalid res_aa_ ( " + string_of( res_aa ) + " )!" );

	utility_exit_with_message( "SHOULD NOT REACH THIS POINT OF THE FUNCTION!!" );

	return CS_RGU_params_; //This is just for prevent compiler warning!

}


} //chemical_shift
} //rna
} //scoring
} //core


