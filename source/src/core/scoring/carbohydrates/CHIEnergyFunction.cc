// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/carbohydrates/CHIEnergyFunction.cc
/// @brief   Method definitions for CHIEnergyFunction.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit Headers
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/database_io.hh>

// Project Header
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Basic Headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <map>
#include <cmath>
#include <string>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.carbohydrates.CHIEnergyFunction" );


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
/// @param    <x>: an angle, in degrees:\n
/// phi (between -180 and 180), if type is ALPHA_LINKS or BETA_LINKS;\n
/// psi (between 0 and 360), if type is _2AX_3EQ_4AX_LINKS, _2EQ_3AX_4EQ_LINKS, ALPHA_1_6_LINKS, or BETA_1_6_LINKS
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
/// @ref      A.K. Nivedha et al. JCTC 2016, 12, 892-901
Energy
CHIEnergyFunction::operator()( CHIEnergyFunctionLinkageType type, core::Angle x ) const {
	return evaluate_function( type, x );
}

/// @details  E'(x) = Sum of -2((x-b)/c)[ae^-((x-b)^2/c)], where the parameters for a, b, and c depend on the linkage
/// type.
/// @param    <x>: an angle, in degrees:\n
/// phi (between -180 and 180), if type is ALPHA_LINKS or BETA_LINKS;\n
/// psi (between 0 and 360), if type is _2AX_3EQ_4AX_LINKS, _2EQ_3AX_4EQ_LINKS, ALPHA_1_6_LINKS, or BETA_1_6_LINKS
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
/// @ref      A.K. Nivedha et al. JCTC 2016, 12, 892-901
Real
CHIEnergyFunction::evaluate_derivative( CHIEnergyFunctionLinkageType type, core::Angle x ) const
{
	if ( type == LINKAGE_NA ) { return 0.0; }

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

	a_.resize( N_LINK_TYPES );
	b_.resize( N_LINK_TYPES );
	c_.resize( N_LINK_TYPES );
	d_.resize( N_LINK_TYPES );

	string const filepath( "scoring/score_functions/carbohydrates/" );
	map< CHIEnergyFunctionLinkageType, string > const filenames = {
		{ ALPHA_LINKS, "CHI_energy_function_for_alpha_linkages.params" },
		{ BETA_LINKS, "CHI_energy_function_for_beta_linkages.params" },
		{ _2AX_3EQ_4AX_LINKS, "CHI_energy_function_for_2ax_3eq_4ax_linkages.params" },
		{ _2EQ_3AX_4EQ_LINKS, "CHI_energy_function_for_2eq_3ax_4eq_linkages.params" },
		{ ALPHA6_LINKS, "CHI_energy_function_for_alpha6_linkages.params" },
		{ BETA6_LINKS, "CHI_energy_function_for_beta6_linkages.params" }
	};

	for ( CHIEnergyFunctionLinkageType type( FIRST_LINK_TYPE ); type <= N_LINK_TYPES; ++type ) {
		map< char, vector1< Real > > params( read_Gaussian_parameters_from_database_file(
				basic::database::full_name( filepath + filenames.at( type ) ) ) );
		a_[ type ] = params[ 'a' ];
		b_[ type ] = params[ 'b' ];
		c_[ type ] = params[ 'c' ];
		d_[ type ] = params[ 'd' ][ 1 ];  // There is only one d value.
	}
}


// Return single CHI Energy Function term, ae^-((x-b)^2/c), for the given type and index.
Energy
CHIEnergyFunction::evaluate_term( CHIEnergyFunctionLinkageType type, uint i, core::Angle x ) const
{
	if ( type == LINKAGE_NA ) { return 0.0; }

	return a_[ type ][ i ] * exp( -( pow( x - b_[ type ][ i ], 2 ) / c_[ type ][ i ] ) );
}

// Sum the individual terms.
Energy
CHIEnergyFunction::evaluate_function( CHIEnergyFunctionLinkageType type, core::Angle x ) const
{
	if ( type == LINKAGE_NA ) { return 0.0; }

	Energy E( 0.0 );
	Size const N( a_[ type ].size() );  // == b_[ type ].size() == c_[ type ].size()

	for ( uint i( 1 ); i <= N; ++i ) {
		E += evaluate_term( type, i, x );
	}

	E += d_[ type ];

	return E;
}


// Sampling Methods ///////////////////////////////////////////////////////////

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
void
CHIEnergyFunction::setup_for_sampling( core::Real step_size ){
	using utility::to_string;

	//Note:
	// Probability from energy: -ln(p)=E -> p = e^-E
	// Used to get rotamer probabilities from energy in PyRosetta Toolkit - need to double check this.

	//Get phi/get psi linkages.
	TR << "Setting up chi sampling" << std::endl;

	utility::vector1< CHIEnergyFunctionLinkageType > phi_linkages( 2 );
	phi_linkages[ 1 ] = ALPHA_LINKS;
	phi_linkages[ 2 ] = BETA_LINKS;

	utility::vector1< CHIEnergyFunctionLinkageType > psi_linkages( 2 );
	psi_linkages[ 1 ] = _2AX_3EQ_4AX_LINKS;
	psi_linkages[ 2 ] = _2EQ_3AX_4EQ_LINKS;

	// Write 2 for loops.  Difference will be -180, 180; 360.
	for ( core::Size i = 1; i <= core::Size(N_LINK_TYPES); ++i ) {
		//TR << "linkage: " << i << std::endl;
		//std::cout << "Angle,Energy,Probability"<<std::endl;

		CHIEnergyFunctionLinkageType linkage_type = static_cast< CHIEnergyFunctionLinkageType >( i );

		CHIDihedralSamplingData sampling_data;

		sampling_data.linkage_type = linkage_type;
		sampling_data.step_size = step_size;

		for ( Angle dih = -180.0; dih <= 180.0; dih+=step_size ) {
			Energy e = evaluate_function(linkage_type, dih);
			Probability prob = std::exp( -e );

			sampling_data.angles.push_back( dih );
			sampling_data.probabilities.push_back( prob );
			//std::cout <<dih << "," << e << "," << prob << std::endl;

		}

		dihedral_sampling_data_[ linkage_type ] = sampling_data;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
CHIEnergyFunction::sampling_data_setup() const {
	if ( dihedral_sampling_data_.empty() ) {
		return false;
	} else {
		return true;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
bool
CHIEnergyFunction::sampling_data_setup( CHIEnergyFunctionLinkageType linkage_type ) const {
	std::map< CHIEnergyFunctionLinkageType, CHIDihedralSamplingData>::const_iterator const_iter;
	const_iter = dihedral_sampling_data_.find( linkage_type );

	if ( const_iter != dihedral_sampling_data_.end() ) {
		return true;
	} else {
		return false;
	}
}

/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
CHIDihedralSamplingData const &
CHIEnergyFunction::get_chi_sampling_data( CHIEnergyFunctionLinkageType linkage_type ) const {
	// Need to go through iterator.  Fuck you C++

	std::map< CHIEnergyFunctionLinkageType, CHIDihedralSamplingData>::const_iterator const_iter;
	const_iter = dihedral_sampling_data_.find( linkage_type );


	if ( const_iter != dihedral_sampling_data_.end() ) {
		return const_iter->second;
	} else {
		std::string m = "CHISampling Data not found for linkage type "+utility::to_string( linkage_type );
		utility_exit_with_message( m );
	}
}


// Helper methods /////////////////////////////////////////////////////////////
// This allows one to use a for loop with LinkageType enum values.
CHIEnergyFunctionLinkageType &
operator++( CHIEnergyFunctionLinkageType & linkage_type )
{
	linkage_type = static_cast< CHIEnergyFunctionLinkageType >( static_cast< int >( linkage_type ) + 1 );
	return linkage_type;
}

}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core
