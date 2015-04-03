// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman
#include <devel/residual_dipolar_coupling/RDC_main.hh>
#include <core/types.hh>
#include <numeric/numeric.functions.hh>

//C++ headers
#include <sstream>
#include <fstream>
#include <map>


using namespace core;

namespace devel {
namespace residual_dipolar_coupling {

void RDC_data::read_RDC_file(
	std::string const & filename
)

{
	std::cout << "Name of input file" << filename << std::endl;
	//	std::map< core::Size, utility::vector1<devel::residual_dipolar_coupling::RDC> > RDC_data_lines;

	std::ifstream infile( filename.c_str() );
	std::string line;

	devel::residual_dipolar_coupling::RDC_data_set my_data_set;

	while( getline( infile, line ) ) {
		std::istringstream line_stream( line );
		std::string atom1, atom2;
		core::Size res1, res2;
		core::Real Jdipolar;
		//		core::Real Reduced_Jdipolar( 0.0001 );
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> Jdipolar;
		if( !line_stream.fail() ) {
			RDC_data_lines[ get_RDC_data_type(atom1,atom2) ].push_back( RDC(get_RDC_data_type(atom1, atom2), res1, res2, Jdipolar ) );
			my_data_set.add_data_line( get_RDC_data_type(atom1, atom2 ), res1, res2, Jdipolar );
		}
	}

	/*
	std::map< core::Size, utility::vector1< devel::residual_dipolar_coupling::RDC > >::iterator it;
	std::vector< devel::residual_dipolar_coupling::RDC >::iterator itt;
	for( it = RDC_data_lines.begin(); it != RDC_data_lines.end(); ++it) {
		for( itt = it->second.begin(); itt != it->second.end(); ++itt ) {
			std::cout << it->first << " " << itt->res1() << " " << itt->res2() << " " << itt->Jdipolar() << " " << itt->Reduced_Jdipolar() << std::endl;
		}
	}
	*/


}

////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////
core::Size RDC_data::get_RDC_data_type(
	std::string const & atom1,
	std::string const & atom2
)
{

	core::Size RDC_type=1;

 	if ( ( atom1 == "N" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "N" ) )
 		RDC_type = 1;
	//****************** FIX THESE LATER !! **************************
	else if ( ( atom1 == "C" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "C" ) )
		RDC_type = 2;
	else if ( ( atom1 == "C" && atom2 == "N" ) || ( atom1 == "N" && atom2 == "C" ) )
		RDC_type = 4;
 	else if ( ( atom1 == "C" && atom2 == "C" ) )
		RDC_type = 5;
	else if ( atom1 == "H" && atom2 == "H" )
		RDC_type = 6;

	//	std::cout << "RDC_type " << RDC_type << std::endl;

	return RDC_type;
}

////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////
core::Real RDC_data::get_invDmax(
	std::string const & atom1,
	std::string const & atom2
)
{

	using namespace numeric;

	core::Size RDC_type( get_RDC_data_type( atom1, atom2 ) );
	core::Real invDcnst( 0.0 );
	core::Real fixed_dist( 0.0 );
	if( RDC_type == 1 ) {
		invDcnst = 0.0000821215;
		fixed_dist = 1.01;
	}
	//************************* FIX THESE LATER !!! ************************
	else if ( RDC_type == 2 ) {
		invDcnst = 0.0000821215;
		fixed_dist = 1.08;
	}
	else if ( RDC_type == 3 ) {
		invDcnst = -0.0000331025;
		fixed_dist = 1.52325877;
	}
	else if ( RDC_type == 4 ) {
		invDcnst = 0.000326533;
		fixed_dist = 1.32874878;
	}
	else if ( RDC_type == 5 ) {
		invDcnst = -0.000131623;
		fixed_dist = 2.032764;
	}
	return invDcnst*cube( fixed_dist );

}


////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////
utility::vector1< devel::residual_dipolar_coupling::RDC  >read_RDC_file(
	std::string const & filename
)

{
	std::cout << "Name of input file" << filename << std::endl;
	//	std::map< core::Size, utility::vector1<devel::residual_dipolar_coupling::RDC> > tmp_RDC_data_lines;
	utility::vector1< devel::residual_dipolar_coupling::RDC > All_RDC_lines;

	std::ifstream infile( filename.c_str() );
	std::string line;


	while( getline( infile, line ) ) {
		std::istringstream line_stream( line );
		std::string atom1, atom2;
		core::Size res1, res2;
		core::Real Jdipolar;
		//		core::Real Reduced_Jdipolar( 0.0001 );
		line_stream >> res1 >> atom1 >> res2 >> atom2 >> Jdipolar;
		if( !line_stream.fail() ) {
			if( res1 == res2 ) {
			core::Size data_type( get_RDC_data_type( atom1, atom2 ) );
			All_RDC_lines.push_back( RDC( data_type, res1, res2, Jdipolar ) );
			} else {
				std::cout << "Skipping this RDC data line. res1 != res2" << std::endl;
			}
		}
	}

	//	All_RDC_sets.push_back( my_data_set );

	utility::vector1< devel::residual_dipolar_coupling::RDC >::iterator it;

	for ( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it ) {
		std::cout << "All RDC lines " << it->type() << " " << it->res1() << " " << it->res2() << " " << it->Jdipolar() << " " << it->fixed_dist() << " " << it->Reduced_Jdipolar() << std::endl;
	}

	return All_RDC_lines;

	/*
	std::map< core::Size, utility::vector1< devel::residual_dipolar_coupling::RDC > >::iterator it;
	std::vector< devel::residual_dipolar_coupling::RDC >::iterator itt;
	for( it = RDC_data_lines.begin(); it != RDC_data_lines.end(); ++it) {
		for( itt = it->second.begin(); itt != it->second.end(); ++itt ) {
			std::cout << it->first << " " << itt->res1() << " " << itt->res2() << " " << itt->Jdipolar() << " " << itt->Reduced_Jdipolar() << std::endl;
		}
	}
	*/
}

////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////

core::Size get_RDC_data_type(
	std::string const & atom1,
	std::string const & atom2
)
{

	core::Size RDC_type( 0 );

 	if ( ( atom1 == "N" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "N" ) )
 		RDC_type = 1;
	//****************** FIX THESE LATER !! **************************
	else if ( ( atom1 == "C" && atom2 == "H" ) || ( atom1 == "H" && atom2 == "C" ) )
		RDC_type = 2;
	else if ( ( atom1 == "C" && atom2 == "N" ) || ( atom1 == "N" && atom2 == "C" ) )
		RDC_type = 4;
 	else if ( ( atom1 == "C" && atom2 == "C" ) )
		RDC_type = 5;
	else if ( atom1 == "H" && atom2 == "H" )
		RDC_type = 6;

	//	std::cout << "RDC_type " << RDC_type << std::endl;

	assert( RDC_type != 0 );

	return RDC_type;
}


////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////


}//ResidualDipolarCoupling
}//devel
