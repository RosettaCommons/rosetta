// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/rna/RNA_SuiteName.cc
/// @brief  RNA suite assignment ported from suitename program (V 0.3.070628)
/// @author  Fang-Chieh Chou

// Unit headers
#include <core/pose/rna/RNA_SuiteName.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>

// Project headers
#include <core/chemical/rna/util.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/fixedsizearray1.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <cmath>

static basic::Tracer TR( "core.pose.rna.RNA_SuiteName" );

namespace core {
namespace pose {
namespace rna {

//Constructor////////////////////////////////
RNA_SuiteName::RNA_SuiteName():
	utility::pointer::ReferenceCount(),
	epsilonmin( 155 ),
	epsilonmax( 310 ),
	delta3min( 55 ),
	delta3max( 110 ),
	delta2min( 120 ),
	delta2max( 175 ),
	gammapmin( 20 ),
	gammapmax( 95 ),
	gammatmin( 140 ),
	gammatmax( 215 ),
	gammammin( 260 ),
	gammammax( 335 ),
	alphamin( 25 ),
	alphamax( 335 ),
	betamin( 50 ),
	betamax( 290 ),
	zetamin( 25 ),
	zetamax( 335 ),
	delta_cutoff( 115 ),
	gamma_cutoff1( 120 ),
	gamma_cutoff2( 240 ),
	suite_undefined_( "__", -1, -1 ),
	dist_pow_( 3 )
{
	init();
}

RNA_SuiteName::~RNA_SuiteName() {}
////////////////////////////////////////////////////////
void RNA_SuiteName::init() {
	using utility::tools::make_vector1;
	using utility::fixedsizearray1;
	//Parameters copied from suitename program//////////////////////
	all_suites_.clear();
	all_suites_.reserve( 54 );
	all_suites_.emplace_back( "1a", 330, fixedsizearray1< Real, 7 >{ { 81.495, 212.250, 288.831, 294.967, 173.990, 53.550, 81.035 } } );
	all_suites_.emplace_back( "1m", 330, fixedsizearray1< Real, 7 >{ { 83.513, 218.120, 291.593, 292.247, 222.300, 58.067, 86.093 } } );
	all_suites_.emplace_back( "1L", 330, fixedsizearray1< Real, 7 >{ { 85.664, 245.014, 268.257, 303.879, 138.164, 61.950, 79.457 } } );
	all_suites_.emplace_back( "&a", 330, fixedsizearray1< Real, 7 >{ { 82.112, 190.682, 264.945, 295.967, 181.839, 51.455, 81.512 } } );
	all_suites_.emplace_back( "7a", 330, fixedsizearray1< Real, 7 >{ { 83.414, 217.400, 222.006, 302.856, 160.719, 49.097, 82.444 } } );
	all_suites_.emplace_back( "3a", 330, fixedsizearray1< Real, 7 >{ { 85.072, 216.324, 173.276, 289.320, 164.132, 45.876, 84.956 } } );
	all_suites_.emplace_back( "9a", 330, fixedsizearray1< Real, 7 >{ { 83.179, 210.347, 121.474, 288.568, 157.268, 49.347, 81.047 } } );
	all_suites_.emplace_back( "1g", 330, fixedsizearray1< Real, 7 >{ { 80.888, 218.636, 290.735, 167.447, 159.565, 51.326, 85.213 } } );
	all_suites_.emplace_back( "7d", 330, fixedsizearray1< Real, 7 >{ { 83.856, 238.750, 256.875, 69.562, 170.200, 52.800, 85.287 } } );
	all_suites_.emplace_back( "3d", 330, fixedsizearray1< Real, 7 >{ { 85.295, 244.085, 203.815, 65.880, 181.130, 54.680, 86.035 } } );
	all_suites_.emplace_back( "5d", 330, fixedsizearray1< Real, 7 >{ { 79.671, 202.471, 063.064, 68.164, 143.450, 49.664, 82.757 } } );
	all_suites_.emplace_back( "1e", 331, fixedsizearray1< Real, 7 >{ { 80.514, 200.545, 280.510, 249.314, 82.662, 167.890, 85.507 } } );
	all_suites_.emplace_back( "1c", 331, fixedsizearray1< Real, 7 >{ { 80.223, 196.591, 291.299, 153.060, 194.379, 179.061, 83.648 } } );
	all_suites_.emplace_back( "1f", 331, fixedsizearray1< Real, 7 >{ { 81.395, 203.030, 294.445, 172.195, 138.540, 175.565, 84.470 } } );
	all_suites_.emplace_back( "5j", 331, fixedsizearray1< Real, 7 >{ { 87.417, 223.558, 80.175, 66.667, 109.150, 176.475, 83.833 } } );
	all_suites_.emplace_back( "5n", 331, fixedsizearray1< Real, 7 >{ { 86.055, 246.502, 100.392, 73.595, 213.752, 183.395, 85.483 } } );
	all_suites_.emplace_back( "1b", 320, fixedsizearray1< Real, 7 >{ { 84.215, 215.014, 288.672, 300.420, 177.476, 58.307, 144.841 } } );
	all_suites_.emplace_back( "1[", 320, fixedsizearray1< Real, 7 >{ { 82.731, 220.463, 288.665, 296.983, 221.654, 54.213, 143.771 } } );
	all_suites_.emplace_back( "3b", 320, fixedsizearray1< Real, 7 >{ { 84.700, 226.400, 168.336, 292.771, 177.629, 48.629, 147.950 } } );
	all_suites_.emplace_back( "1z", 320, fixedsizearray1< Real, 7 >{ { 83.358, 206.042, 277.567, 195.700, 161.600, 50.750, 145.258 } } );
	all_suites_.emplace_back( "5z", 320, fixedsizearray1< Real, 7 >{ { 82.614, 206.440, 52.524, 163.669, 148.421, 50.176, 147.590 } } );
	all_suites_.emplace_back( "7p", 320, fixedsizearray1< Real, 7 >{ { 84.285, 236.600, 220.400, 68.300, 200.122, 53.693, 145.730 } } );
	all_suites_.emplace_back( "5p", 320, fixedsizearray1< Real, 7 >{ { 84.457, 213.286, 69.086, 75.500, 156.671, 57.486, 147.686 } } );
	all_suites_.emplace_back( "1t", 321, fixedsizearray1< Real, 7 >{ { 81.200, 199.243, 288.986, 180.286, 194.743, 178.200, 147.386 } } );
	all_suites_.emplace_back( "5q", 321, fixedsizearray1< Real, 7 >{ { 82.133, 204.933, 69.483, 063.417, 115.233, 176.283, 145.733 } } );
	all_suites_.emplace_back( "1o", 322, fixedsizearray1< Real, 7 >{ { 83.977, 216.508, 287.192, 297.254, 225.154, 293.738, 150.677 } } );
	all_suites_.emplace_back( "7r", 322, fixedsizearray1< Real, 7 >{ { 84.606, 232.856, 248.125, 63.269, 181.975, 295.744, 149.744 } } );
	all_suites_.emplace_back( "5r", 322, fixedsizearray1< Real, 7 >{ { 83.000, 196.900, 65.350, 60.150, 138.425, 292.550, 154.275 } } );
	all_suites_.emplace_back( "2a", 230, fixedsizearray1< Real, 7 >{ { 145.399, 260.339, 288.756, 288.444, 192.733, 53.097, 84.067 } } );
	all_suites_.emplace_back( "4a", 230, fixedsizearray1< Real, 7 >{ { 146.275, 259.783, 169.958, 298.450, 169.583, 50.908, 83.967 } } );
	all_suites_.emplace_back( "0a", 230, fixedsizearray1< Real, 7 >{ { 149.286, 223.159, 139.421, 284.559, 158.107, 47.900, 84.424 } } );
	all_suites_.emplace_back( "#a", 230, fixedsizearray1< Real, 7 >{ { 148.006, 191.944, 146.231, 289.288, 150.781, 42.419, 84.956 } } );
	all_suites_.emplace_back( "4g", 230, fixedsizearray1< Real, 7 >{ { 148.028, 256.922, 165.194, 204.961, 165.194, 49.383, 82.983 } } );
	all_suites_.emplace_back( "6g", 230, fixedsizearray1< Real, 7 >{ { 145.337, 262.869, 79.588, 203.863, 189.688, 58.000, 84.900 } } );
	all_suites_.emplace_back( "8d", 230, fixedsizearray1< Real, 7 >{ { 148.992, 270.596, 240.892, 62.225, 176.271, 53.600, 87.262 } } );
	all_suites_.emplace_back( "4d", 230, fixedsizearray1< Real, 7 >{ { 149.822, 249.956, 187.678, 80.433, 198.133, 61.000, 89.378 } } );
	all_suites_.emplace_back( "6d", 230, fixedsizearray1< Real, 7 >{ { 146.922, 241.222, 88.894, 59.344, 160.683, 052.333, 83.417 } } );
	all_suites_.emplace_back( "2g", 230, fixedsizearray1< Real, 7 >{ { 141.900, 258.383, 286.517, 178.267, 165.217, 48.350, 84.783 } } );
	all_suites_.emplace_back( "2h", 231, fixedsizearray1< Real, 7 >{ { 147.782, 260.712, 290.424, 296.200, 177.282, 175.594, 86.565 } } );
	all_suites_.emplace_back( "4n", 231, fixedsizearray1< Real, 7 >{ { 143.722, 227.256, 203.789, 73.856, 216.733, 194.444, 80.911 } } );
	all_suites_.emplace_back( "0i", 231, fixedsizearray1< Real, 7 >{ { 148.717, 274.683, 100.283, 80.600, 248.133, 181.817, 82.600 } } );
	all_suites_.emplace_back( "6n", 231, fixedsizearray1< Real, 7 >{ { 150.311, 268.383, 84.972, 63.811, 191.483, 176.644, 85.600 } } );
	all_suites_.emplace_back( "6j", 231, fixedsizearray1< Real, 7 >{ { 141.633, 244.100, 66.056, 71.667, 122.167, 182.200, 83.622 } } );
	all_suites_.emplace_back( "0k", 232, fixedsizearray1< Real, 7 >{ { 149.070, 249.780, 111.520, 278.370, 207.780, 287.820, 86.650 } } );
	all_suites_.emplace_back( "2[", 220, fixedsizearray1< Real, 7 >{ { 146.383, 259.402, 291.275, 291.982, 210.048, 54.412, 147.760 } } );
	all_suites_.emplace_back( "4b", 220, fixedsizearray1< Real, 7 >{ { 145.256, 244.622, 162.822, 294.159, 171.630, 45.900, 145.804 } } );
	all_suites_.emplace_back( "0b", 220, fixedsizearray1< Real, 7 >{ { 147.593, 248.421, 112.086, 274.943, 164.764, 56.843, 146.264 } } );
	all_suites_.emplace_back( "4p", 220, fixedsizearray1< Real, 7 >{ { 150.077, 260.246, 213.785, 71.900, 207.638, 56.715, 148.131 } } );
	all_suites_.emplace_back( "6p", 220, fixedsizearray1< Real, 7 >{ { 146.415, 257.831, 89.597, 67.923, 173.051, 55.513, 147.623 } } );
	all_suites_.emplace_back( "2z", 220, fixedsizearray1< Real, 7 >{ { 142.900, 236.550, 268.800, 180.783, 185.133, 54.467, 143.350 } } );
	all_suites_.emplace_back( "4s", 221, fixedsizearray1< Real, 7 >{ { 149.863, 247.562, 170.488, 277.938, 84.425, 176.413, 148.087 } } );
	all_suites_.emplace_back( "2u", 221, fixedsizearray1< Real, 7 >{ { 143.940, 258.200, 298.240, 279.640, 183.680, 183.080, 145.120 } } );
	all_suites_.emplace_back( "2o", 222, fixedsizearray1< Real, 7 >{ { 147.342, 256.475, 295.508, 287.408, 194.525, 293.725, 150.458 } } );

	utility::vector1<std::string> _dominant_suites  = make_vector1( "1a", "1c", "1b", "0a", "6n" );
	utility::vector1<std::string> _satellite_suites = make_vector1( "1m", "1L", "&a", "1f", "1[", "4a", "#a", "0i", "6j" );

	regular_half_width_ = fixedsizearray1< Size, 7 >{ { 28, 60, 55, 50, 70, 35, 28 } };
	dominant_suites_ = _dominant_suites;
	satellite_suites_ = _satellite_suites;
	half_width_sat_ = make_vector1(
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 32, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 18, 55, 50, 18, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 20, 20, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 47, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 34, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 40, 40, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 26, 26, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 70, 35, 28 } } );

	half_width_dom_ = make_vector1(
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 64, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 70, 55, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 60, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 65, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 56, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 50, 50, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 36, 36, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 70, 35, 28 } },
		fixedsizearray1< Size, 7 > { { 28, 60, 55, 50, 70, 35, 28 } } );
}

RNA_SuiteInfo
RNA_SuiteName::closest_suite( Pose const & pose, Size const res ) const {
	using namespace chemical::rna;
	if ( is_rna_chainbreak( pose, res - 1 ) ) return closest_by_dist4( pose, res );
	utility::fixedsizearray1< Real, 7 > torsions{ {
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( DELTA ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( EPSILON ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( ZETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( ALPHA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( BETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( GAMMA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( DELTA ) ) } };
	return closest_suite( torsions );
}

RNA_SuiteInfo
RNA_SuiteName::closest_by_dist4( Pose const & pose, Size const res ) const {
	using namespace chemical::rna;
	utility::fixedsizearray1<Real, 7> torsions{ {
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( DELTA ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( EPSILON ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( ZETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( ALPHA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( BETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( GAMMA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( DELTA ) ) } };
	return closest_by_dist4( torsions );
}

RNA_SuiteInfo
RNA_SuiteName::closest_suite(
	utility::fixedsizearray1< Real, 7 > const & res_dihedrals
) const {
	RNA_SuiteInfo closest;
	Real distance = -1;
	for ( auto const & suite : all_suites_ ) {
		Real const this_dist = distance_7d( suite.torsion, res_dihedrals, regular_half_width_ );
		if ( distance < 0 || this_dist < distance ) {
			distance = this_dist;
			closest = suite;
		}
	}
	TR << "Closest: resid torsions " << res_dihedrals << std::endl;
	TR << "Closest: suite torsions " << closest.torsion << std::endl;
	return closest;
}

RNA_SuiteInfo
RNA_SuiteName::closest_by_dist4(
	utility::fixedsizearray1< Real, 7 > const & res_dihedrals
) const {
	RNA_SuiteInfo closest;
	Real distance = -1;
	for ( auto const & suite : all_suites_ ) {
		Real const this_dist = distance_4d( suite.torsion, res_dihedrals, regular_half_width_ );
		if ( distance < 0 || this_dist < distance ) {
			distance = this_dist;
			closest = suite;
		}
	}
	TR << "Closest: resid torsions " << res_dihedrals << std::endl;
	TR << "Closest: suite torsions " << closest.torsion << std::endl;
	return closest;
}

/////////////////////////////////////////////////////////
// Called by RNA_SuitePotential which uses this RNA_SuiteName
//  class to calculate suiteness. The user may specify a suiteness_bonus file
//  with alternative rotamer centers -- or new ones.
// Note use of string tags. May be slow, but allows us to preserve Richardson
//  specification of a few 'dominant' clusters and satellite clusters.
void
RNA_SuiteName::update_centers( utility::vector1< utility::fixedsizearray1< Real, 7 > > const & centers,
	utility::vector1< std::string > const & tags ){
	utility::vector1< RNA_SuiteInfo > new_suites;
	for ( Size n = 1; n <= tags.size(); n++ ) {
		bool found_match = false;
		for ( auto & suite : all_suites_ ) {
			if ( suite.name == tags[ n ] ) {
				found_match = true;
				suite.torsion = centers[ n ];
				new_suites.push_back( suite );
			}
		}

		if ( !found_match ) {
			std::string const & tag = tags[n];
			runtime_assert( tag.size() > 0 );
			new_suites.push_back( RNA_SuiteInfo( tag, get_classifier( centers[n] ), centers[n] ) );
		}
	}
	all_suites_ = new_suites;
}

/////////////////////////////////////////////////////////
RNA_SuiteInfo
RNA_SuiteName::name2suite( std::string const & name ) const{
	for ( auto const & suite : all_suites_ ) {
		if ( suite.name == name ) {
			return suite;
		}
	}
	utility_exit_with_message( "Invalid suitename!" );
	return all_suites_[1];
}

class principal_angle_degrees {
public:
	Real operator()( Real r ) {
		return numeric::nonnegative_principal_angle_degrees( r );
	}
};

///////////////////////////////////
//Suite assignment codes
Size
RNA_SuiteName::get_classifier( utility::fixedsizearray1< Real,7 > const & torsions,
	bool & is_outlier ) const {

	std::for_each( torsions.begin(), torsions.end(), principal_angle_degrees() );

	Size classifier = 0;
	//Delta1
	if ( torsions[1] >= delta3min && torsions[1] <= delta3max ) {
		classifier += 300;
	} else if ( torsions[1] >= delta2min && torsions[1] <= delta2max ) {
		classifier += 200;
	} else {
		is_outlier = true;
		if ( torsions[1] <= delta_cutoff ) {
			classifier += 300;
		} else {
			classifier += 200;
		}
	}
	//Delta2
	if ( torsions[7] >= delta3min && torsions[7] <= delta3max ) {
		classifier += 30;
	} else if ( torsions[7] >= delta2min && torsions[7] <= delta2max ) {
		classifier += 20;
	} else {
		is_outlier = true;
		if ( torsions[7] <= delta_cutoff ) {
			classifier += 30;
		} else {
			classifier += 20;
		}
	}
	//Gamma
	if ( torsions[6] >= gammapmin && torsions[6] <= gammapmax ) {
		;
	} else if ( torsions[6] >= gammatmin && torsions[6] <= gammatmax ) {
		classifier += 1;
	} else if ( torsions[6] >= gammammin && torsions[6] <= gammammax ) {
		classifier += 2;
	} else {
		is_outlier = true;
		if ( torsions[7] <= gamma_cutoff1 ) {
			;
		} else if ( torsions[7] > gamma_cutoff2 ) {
			classifier += 2;
		} else {
			classifier += 1;
		}
	}
	return classifier;
}

////////////////////////////////////////////////////////////
Size
RNA_SuiteName::get_classifier( utility::fixedsizearray1< Real,7 > const & torsions ) const {
	bool is_outlier; // dummy
	return get_classifier( torsions, is_outlier );
}

//Distance computation///////////////////////////////
Real RNA_SuiteName::distance_4d(
	utility::fixedsizearray1<Real,7> const & torsion1,
	utility::fixedsizearray1<Real,7> const & torsion2,
	utility::fixedsizearray1<Size, 7> const & half_width
) const {
	// By suitename default, distance with power of 3 is used.
	Real sum = 0;
	for ( Size i = 2; i <= 5; ++i ) {
		Real diff = numeric::principal_angle_degrees( torsion1[i] - torsion2[i] );
		diff = diff / static_cast<Real>( half_width[i] );
		if ( diff < 0 ) diff = -diff;
		sum += diff * diff * diff;
	}

	return pow(sum, 1.0 / dist_pow_ );
}

////////////////////////////////////////////
Real RNA_SuiteName::distance_7d(
	utility::fixedsizearray1<Real,7> const & torsion1,
	utility::fixedsizearray1<Real,7> const & torsion2,
	utility::fixedsizearray1<Size, 7> const & half_width
) const {
	utility::fixedsizearray1<Real,7> deriv;
	return distance_7d( torsion1, torsion2, half_width, deriv );
}

////////////////////////////////////////////
Real RNA_SuiteName::distance_7d(
	utility::fixedsizearray1<Real,7> const & torsion1,
	utility::fixedsizearray1<Real,7> const & torsion2,
	utility::fixedsizearray1<Size, 7> const & half_width,
	utility::fixedsizearray1<Real,7> & deriv
) const {
	//By suitename default, distance with power of 3 is used.
	Real sum = 0;
	for ( Size i = 1; i <= 7; ++i ) {
		Real diff = numeric::principal_angle_degrees( torsion1[i] - torsion2[i] );
		diff = diff / static_cast<Real>( half_width[i] );
		int sign = ( diff < 0 ) ? -1 : 1;
		if ( sign < 0 ) diff = -diff;
		sum += pow( diff, dist_pow_ );

		deriv[ i ] = sign * pow( diff, (dist_pow_ - 1) ) / static_cast<Real>(half_width[i]) ;
	}
	for ( Size i = 1; i <= deriv.size(); ++i ) deriv[ i ] *=  pow( sum, (1.0 / dist_pow_ ) - 1 );

	return pow( sum, 1.0 / dist_pow_ );
}

//////////////////////////////////////////////
bool RNA_SuiteName::is_in_between(
	utility::fixedsizearray1<Real,7> const & target,
	utility::fixedsizearray1<Real,7> const & dominant,
	utility::fixedsizearray1<Real,7> const & satellite
) const {
	Real inner_product_dom( 0 ), inner_product_sat( 0 );
	for ( Size i = 2; i <= 5; ++i ) {
		inner_product_dom += ( target[i] - dominant[i]) * (satellite[i] - dominant[i] );
		inner_product_sat += ( target[i] - satellite[i]) * (dominant[i] - satellite[i] );
	}
	if ( inner_product_dom < 0 ) return false;
	if ( inner_product_sat < 0 ) return false;
	return true;
}

////////////////////////////////
Real
RNA_SuiteName::get_suiteness( Real const & dist_7d ) const {
	runtime_assert( dist_7d <= 1.0 );
	return ( cos( numeric::constants::r::pi * dist_7d ) + 1.0 ) * 0.5;
}

////////////////////////////////
Real
RNA_SuiteName::get_suiteness_derivative( Real const & dist_7d ) const {
	runtime_assert( dist_7d <= 1.0 );
	return ( -0.5 * numeric::constants::r::pi * sin( numeric::constants::r::pi * dist_7d ) );
}

////////////////////////////////
void
RNA_SuiteName::fill_suiteness_derivative_7d (
	utility::fixedsizearray1< Real, 7 > & deriv,
	utility::fixedsizearray1< Real,7 > const & torsions,
	utility::fixedsizearray1< Real,7 > const & torsions_center,
	utility::fixedsizearray1< Size,7 > const & half_width ) const {

	Distance const dist7 = distance_7d( torsions, torsions_center, half_width, deriv );

	// chain rule.
	Real const suite_derivative = get_suiteness_derivative( dist7 );
	for ( Size n = 1; n <= torsions.size(); n++ ) deriv[ n ] *= suite_derivative;
}

//Suite assign/////////////////////////
RNA_SuiteAssignment
RNA_SuiteName::assign( Pose const & pose, Size const res ) const {
	using namespace chemical::rna;
	if ( res == 1 ) return suite_undefined_;
	if ( res > pose.total_residue() ) return suite_undefined_;
	if ( !pose.residue( res - 1 ).is_RNA()  ) return suite_undefined_;
	if ( !pose.residue( res ).is_RNA()  )     return suite_undefined_;
	if ( is_rna_chainbreak(pose, res - 1) )   return suite_undefined_;
	if ( is_rna_chainbreak(pose, res - 1) ) return suite_undefined_;
	utility::fixedsizearray1< Real, 7 > torsions{ {
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( DELTA ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( EPSILON ) ),
		numeric::principal_angle_degrees( pose.residue( res - 1 ).mainchain_torsion( ZETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( ALPHA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( BETA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( GAMMA ) ),
		numeric::principal_angle_degrees( pose.residue( res ).mainchain_torsion( DELTA ) ) } };
	return assign( torsions );
}

////////////////////////////
RNA_SuiteAssignment
RNA_SuiteName::assign( utility::fixedsizearray1<Real,7> const & torsions_in ) const {
	utility::fixedsizearray1< Real, 7 > deriv; // dummy. empty vector means: don't compute derivs.
	return assign( torsions_in, deriv );
}

////////////////////////////

RNA_SuiteAssignment
RNA_SuiteName::assign( utility::fixedsizearray1<Real,7> const & torsions_in,
	utility::fixedsizearray1<Real,7> & deriv ) const {
	using namespace chemical::rna;
	debug_assert( torsions_in.size() == 7 );

	utility::fixedsizearray1<Real,7> torsions = torsions_in;
	for ( Real & torsion : torsions ) {
		torsion = std::fmod( torsion, 360.0 );
		if ( torsion < 0 ) torsion += 360;
	}
	for ( Real & d : deriv ) d = 0;

	bool is_outlier( false );

	//Fast check for outlier torsions
	if ( torsions[2] < epsilonmin || torsions[2] > epsilonmax ) is_outlier = true;
	if ( torsions[3] < zetamin || torsions[3] > zetamax ) is_outlier = true;
	if ( torsions[4] < alphamin || torsions[4] > alphamax ) is_outlier = true;
	if ( torsions[5] < betamin || torsions[5] > betamax ) is_outlier = true;

	//Classify the torsion using delta and gamma
	Size const classifier = get_classifier( torsions, is_outlier );

	// Special case: 332 has no valid suite
	if ( classifier == 332 ) {
		RNA_SuiteAssignment const outlier( "!!", 0, 5 );
		return outlier;
	}

	Size best_index( 0 ), dom_index( 0 );
	Real best_dist( -1 ), dom_dist( 999 );

	for ( Size i = 1, e = all_suites_.size(); i <= e; ++i ) {
		auto const & suite = all_suites_[ i ];
		if ( suite.classifier != classifier ) continue;
		Real const dist = distance_4d( torsions, suite.torsion,
			regular_half_width_ );
		if ( dist <= 1 && dominant_suites_.has_value( suite.name ) ) {
			dom_index = i;
			dom_dist = dist;
		} else if ( dist < best_dist || best_dist < 0 ) {
			best_dist = dist;
			best_index = i;
		}
	}

	//std::cout << best_index << ' ' << dom_index << ' ' << classifier <<std::endl;
	if ( best_dist > 1 && dom_dist <= 1 ) {
		best_index = dom_index;
	} else {
		if ( dom_index != 0 ) {
			Size const find_index = satellite_suites_.index(
				all_suites_[best_index].name );
			if ( satellite_suites_.has_value(all_suites_[best_index].name ) &&
					is_in_between( torsions, all_suites_[dom_index].torsion,
					all_suites_[best_index].torsion ) ) {
				Real const satellite_dist = distance_4d( torsions,
					all_suites_[best_index].torsion, half_width_sat_[find_index] );
				Real const dominant_dist = distance_4d( torsions,
					all_suites_[dom_index].torsion, half_width_dom_[find_index] );
				if ( satellite_dist > dominant_dist ) {
					best_index = dom_index;
				}
			} else if ( best_dist > dom_dist ) {
				best_index = dom_index;
			}
		}
	}

	Real const dist_7d = distance_7d( torsions, all_suites_[best_index].torsion, regular_half_width_ );
	if ( dist_7d > 1 || is_outlier ) {
		RNA_SuiteAssignment const outlier( "!!", 0, dist_7d );
		return outlier;
	}

	Real const suiteness = get_suiteness( dist_7d );

	fill_suiteness_derivative_7d( deriv, torsions, all_suites_[best_index].torsion, regular_half_width_ );

	std::string const suitename = all_suites_[best_index].name;
	RNA_SuiteAssignment const best_suite ( suitename, suiteness, dist_7d );
	return best_suite;
}


}
}
}
