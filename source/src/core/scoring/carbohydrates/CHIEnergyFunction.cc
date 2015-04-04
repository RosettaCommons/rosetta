// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/carbohydrates/CHIEnergyFunction.cc
/// @brief   Method definitions for CHIEnergyFunction.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/database_io.hh>

// Project Header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// Basic header
#include <basic/database/open.hh>

// C++ headers
#include <map>
#include <cmath>
#include <string>


namespace core {
namespace scoring {
namespace carbohydrates {

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Default constructor
/// @details  This class is only intended to be instantiated by the ScoringManager.
CHIEnergyFunction::CHIEnergyFunction()
{
	init();
}

CHIEnergyFunction::~CHIEnergyFunction()
{}


// Other Public Methods ///////////////////////////////////////////////////////
/// @details  E(x) = d + Sum of ae^-((x-b)^2/c), where the parameters for a, b, c, and d depend on the linkage type.
/// @param    <x>: an angle, in degrees\n
/// phi (between -180 and 180), if type is ALPHA_LINKS or BETA_LINKS;\n
/// psi (between 0 and 360), if type is _2AX_3EQ_4AX_LINKS or _2EQ_3AX_4EQ_LINKS
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
Energy
CHIEnergyFunction::operator()( LinkageType type, core::Angle x ) const {
	return evaluate_function( type, x );
}

/// @details  E'(x) = Sum of -2((x-b)/c)[ae^-((x-b)^2/c)], where the parameters for a, b, and c depend on the linkage
/// type.
/// @param    <x>: an angle, in degrees, between 0 and 360:\n
/// phi, if type is ALPHA_LINKS or BETA_LINKS; psi, if type is _2AX_3EQ_4AX_LINKS or _2EQ_3AX_4EQ_LINKS
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
Real
CHIEnergyFunction::evaluate_derivative( LinkageType type, core::Angle x ) const
{
	Real derivative( 0.0 );
	Size const N( a_[ type ].size() );  // == b_[ type ].size() == c_[ type ].size()

	for ( uint i( 1 ); i <= N; ++i ) {
		derivative += -2 * ( ( x - b_[ type ][ i ] ) / c_[ type ][ i ] ) * evaluate_term( type, i, x );
	}

	return derivative;
}


// Private methods ////////////////////////////////////////////////////////////
void
CHIEnergyFunction::init()
{
	using namespace std;
	using namespace utility;

	//std::map< char, utility::vector1< Real > > read_Gaussian_parameters_from_database_file( std::string const & filename );
	a_.resize( N_LINK_TYPES );
	b_.resize( N_LINK_TYPES );
	c_.resize( N_LINK_TYPES );
	d_.resize( N_LINK_TYPES );

	for ( LinkageType type( FIRST_LINK_TYPE ); type <= N_LINK_TYPES; ++type ) {
		string const filepath( "scoring/score_functions/carbohydrates/" );
		string filename;
		switch ( type ) {
			case ALPHA_LINKS:
				filename = "CHI_energy_function_for_alpha_linkages.params";
				break;
			case BETA_LINKS:
				filename = "CHI_energy_function_for_beta_linkages.params";
				break;
			case _2AX_3EQ_4AX_LINKS:
				filename = "CHI_energy_function_for_2ax_3eq_4ax_linkages.params";
				break;
			case _2EQ_3AX_4EQ_LINKS:
				filename = "CHI_energy_function_for_2eq_3ax_4eq_linkages.params";
				break;
		}

		map< char, vector1< Real > > params( read_Gaussian_parameters_from_database_file(
				basic::database::full_name( filepath + filename ) ) );
		a_[ type ] = params[ 'a' ];
		b_[ type ] = params[ 'b' ];
		c_[ type ] = params[ 'c' ];
		d_[ type ] = params[ 'd' ][ 1 ];  // There is only one d value.
	}
}


// Return single CHI energy function term, ae^-((x-b)^2/c), for the given type and index.
Energy
CHIEnergyFunction::evaluate_term( LinkageType type, uint i, core::Angle x ) const
{
	return a_[ type ][ i ] * exp( -( pow( x - b_[ type ][ i ], 2 ) / c_[ type ][ i ] ) );
}

// Sum the individual terms.
Energy
CHIEnergyFunction::evaluate_function( LinkageType type, core::Angle x ) const
{
	Energy E( 0.0 );
	Size const N( a_[ type ].size() );  // == b_[ type ].size() == c_[ type ].size()

	for ( uint i( 1 ); i <= N; ++i ) {
		E += evaluate_term( type, i, x );
	}

	E += d_[ type ];

	return E;
}


// Helper methods /////////////////////////////////////////////////////////////
// This allows one to use a for loop with LinkageType enum values.
LinkageType &
operator++( LinkageType & type )
{
	type = static_cast< LinkageType >( static_cast< int >( type ) + 1 );
	return type;
}

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core
