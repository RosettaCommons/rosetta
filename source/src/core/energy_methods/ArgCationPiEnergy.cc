// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/energy_methods/ArgCationPiEnergy.cc
/// @brief  Cation pi term that specializes in bringing Arginine and rings together
/// @details Currently designed for canonical amino acids but easily extended beyond.
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <core/energy_methods/ArgCationPiEnergy.hh>
#include <core/energy_methods/ArgCationPiEnergyCreator.hh>


// Package headers

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic Headers
#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


static basic::Tracer TR("core.energy_methods.ArgCationPiEnergy");

namespace core {
namespace scoring {
namespace methods {



scoring::methods::EnergyMethodOP
ArgCationPiEnergyCreator::create_energy_method(
	scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< ArgCationPiEnergy >( options );
}

scoring::ScoreTypes
ArgCationPiEnergyCreator::score_types_for_method() const {
	scoring::ScoreTypes sts;
	sts.push_back( arg_cation_pi );
	return sts;
}

ArgCationPiEnergy::ArgCationPiEnergy( core::scoring::methods::EnergyMethodOptions const & options ) :
	scoring::methods::ContextIndependentTwoBodyEnergy( utility::pointer::make_shared< ArgCationPiEnergyCreator >() ),
	his_can_be_pi_( options.arg_cation_pi_his_can_be_pi() )
{}

scoring::methods::EnergyMethodOP
ArgCationPiEnergy::clone() const {
	return utility::pointer::make_shared< ArgCationPiEnergy >( *this );
}

bool
ArgCationPiEnergy::is_arg(
	core::chemical::ResidueType const & restype
) const {
	using namespace core::chemical;
	AA const aa( restype.aa() );
	return aa == aa_arg || aa == aa_dar || aa == aa_b3r;
}

// Someday when this is generalized to all rings, these functions can go away...

bool
ArgCationPiEnergy::is_phe(
	core::chemical::ResidueType const & restype
) const {
	using namespace core::chemical;
	AA const aa( restype.aa() );
	return aa == aa_phe || aa == aa_dph || aa == aa_b3f;
}

bool
ArgCationPiEnergy::is_his(
	core::chemical::ResidueType const & restype
) const {
	using namespace core::chemical;
	AA const aa( restype.aa() );
	return aa == aa_his || aa == aa_dhi || aa == aa_b3h;
}

bool
ArgCationPiEnergy::is_trp(
	core::chemical::ResidueType const & restype
) const {
	using namespace core::chemical;
	AA const aa( restype.aa() );
	return aa == aa_trp || aa == aa_dtr || aa == aa_b3w;
}

bool
ArgCationPiEnergy::is_tyr(
	core::chemical::ResidueType const & restype
) const {
	using namespace core::chemical;
	AA const aa( restype.aa() );
	return aa == aa_tyr || aa == aa_dty || aa == aa_b3y;
}

bool
ArgCationPiEnergy::is_pi(
	core::chemical::ResidueType const & restype
) const {
	return is_phe( restype ) || is_trp( restype ) || is_tyr( restype ) || ( his_can_be_pi_ && is_his(restype) );
}

bool
ArgCationPiEnergy::valid_res_pair(
	core::chemical::ResidueType const & restype1,
	core::chemical::ResidueType const & restype2
) const {
	if ( is_arg( restype1 ) && is_pi( restype2 ) ) return true;
	if ( is_arg( restype2 ) && is_pi( restype1 ) ) return true;
	return false;
}

std::pair< std::pair< numeric::xyzVector<Real>, numeric::xyzVector<Real> >, utility::vector1<std::string> >
ArgCationPiEnergy::get_ring( conformation::Residue const & pi, bool trp_little_ring ) const {

	utility::vector1<std::string> atom_names;

	if ( is_trp( pi.type() ) ) {
		if ( trp_little_ring ) {
			atom_names = utility::vector1<std::string> { "CG", "CD1", "CD2", "NE1", "CE2" };
		} else {
			atom_names = utility::vector1<std::string> { "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2" };
		}
	}
	if ( is_tyr( pi.type() ) ) {
		atom_names = utility::vector1<std::string> { "CG", "CD1", "CD2", "CE1", "CE2", "CZ" };
	}
	if ( is_phe( pi.type() ) ) {
		atom_names = utility::vector1<std::string> { "CG", "CD1", "CD2", "CE1", "CE2", "CZ" };
	}
	if ( is_his( pi.type() ) ) {
		atom_names = utility::vector1<std::string> { "CG", "ND1", "CD2", "CE1", "NE2" };
	}

	runtime_assert( atom_names.size() > 0 );

	utility::vector1<numeric::xyzVector<Real>> pts;
	numeric::xyzVector<Real> center( 0, 0, 0 );
	for ( std::string const & name : atom_names ) {
		pts.push_back( pi.xyz( name ) );
		center += pts.back();
	}

	center /= pts.size();

	numeric::xyzVector<Real> normal = (pts[1] - center).cross(pts[2] - center).normalized();

	return {{ center, normal }, atom_names};
}


std::pair< std::pair< numeric::xyzVector<Real>, numeric::xyzVector<Real> >, utility::vector1<std::string> >
ArgCationPiEnergy::get_ring_params( conformation::Residue const & pi, conformation::Residue const & arg ) const {

	auto params = get_ring( pi, true );

	if ( is_trp( pi.type() ) ) {
		auto big_ring_params =  get_ring( pi, false );

		numeric::xyzVector<Real> arg_cz = arg.xyz("CZ");

		if ( params.first.first.distance(arg_cz) > big_ring_params.first.first.distance(arg_cz) ) {
			params = big_ring_params;
		}
	}

	return params;
}

Real
point_line_distance(
	numeric::xyzVector<Real> const & point,
	numeric::xyzVector<Real> const & line_base,
	numeric::xyzVector<Real> const & line_ray
) {
	numeric::xyzVector<Real> line_other = line_ray + line_base;

	//http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	return ( point - line_base ).cross( point - line_other ).norm() / line_ray.norm();
}

void
ArgCationPiEnergy::residue_pair_energy(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( ! valid_res_pair( res1.type(), res2.type() ) ) return;

	bool arg_is_1 = is_arg( res1.type() );

	conformation::Residue const & arg = arg_is_1 ? res1 : res2;
	conformation::Residue const & pi = arg_is_1 ? res2 : res1;

	auto params = get_ring_params( pi, arg );
	auto pi_center_normal = params.first;

	numeric::xyzVector<Real> arg_cz = arg.xyz("CZ");

	numeric::xyzVector<Real> arg_normal = ( arg.xyz("NH1") - arg_cz).cross( arg.xyz("NH2") - arg_cz ).normalized();

	Real cylinder_offset = point_line_distance( pi_center_normal.first, arg_cz, arg_normal );

	// You really just have to plot these functions to see what is happening here
	Real cylinder_score = std::max<Real>(0, std::cos( cylinder_offset / 2 ));
	if ( cylinder_offset > 4 ) {
		cylinder_score = 0;
	}

	Real plane_angle = std::acos( std::abs( arg_normal.dot( pi_center_normal.second ) ) );
	Real plane_score = 1;
	if ( plane_angle > numeric::conversions::radians( 10.0 ) ) {
		if ( plane_angle > numeric::conversions::radians( 70.0 ) ) {
			plane_score = 0;
		} else {
			plane_score = std::max<Real>( 0, ( std::cos( 3 * ( plane_angle - numeric::conversions::radians(10.0) ) ) + 1 ) / 2 );
			plane_score *= plane_score;
		}
	}

	Real distance = arg_cz.distance( pi_center_normal.first );
	Real distance_score = 1;
	if ( distance > 3.7 && distance < 6 ) {
		distance_score = std::max<Real>( 0, std::cos( ( distance - 3.7 ) / 1.4 ) );
	}

	Real score = -1 * cylinder_score * plane_score * distance_score;

	emap[ arg_cation_pi ] += score;
}


bool
ArgCationPiEnergy::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool
) const {
	return valid_res_pair( res1.type(), res2.type() );
}


void
ArgCationPiEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const {

	// For now we're only going to do the distance derivative

	if ( ! valid_res_pair( res1.type(), res2.type() ) ) return;

	bool arg_is_1 = is_arg( res1.type() );

	conformation::Residue const & arg = arg_is_1 ? res1 : res2;
	conformation::Residue const & pi = arg_is_1 ? res2 : res1;

	auto params = get_ring_params( pi, arg );
	auto pi_center_normal = params.first;
	auto ring_atoms = params.second;

	numeric::xyzVector<Real> arg_cz = arg.xyz("CZ");

	numeric::xyzVector<Real> arg_normal = ( arg.xyz("NH1") - arg_cz).cross( arg.xyz("NH2") - arg_cz ).normalized();

	Real cylinder_offset = point_line_distance( pi_center_normal.first, arg_cz, arg_normal );


	// You really just have to plot these functions to see what is happening here
	Real cylinder_score = std::max<Real>(0, std::cos( cylinder_offset / 2 ));
	if ( cylinder_offset > 4 ) {
		cylinder_score = 0;
	}

	Real plane_angle = std::acos( std::abs( arg_normal.dot( pi_center_normal.second ) ) );
	Real plane_score = 1;
	if ( plane_angle > numeric::conversions::radians( 10.0 ) ) {
		// std::cout << plane_angle << " vs " <<  numeric::conversions::radians( 70.0 ) << std::endl;
		if ( plane_angle > numeric::conversions::radians( 70.0 ) ) {
			// std::cout << "is zero" << std::endl;
			plane_score = 0;
		} else {
			plane_score = std::max<Real>( 0, ( std::cos( 3 * ( plane_angle - numeric::conversions::radians(10.0) ) ) + 1 ) / 2 );
			// std::cout << ( std::cos( 3 * ( plane_angle - numeric::conversions::radians(10.0) ) ) + 1 ) / 2;
			plane_score *= plane_score;
		}
	}

	Real distance = arg_cz.distance( pi_center_normal.first );
	Real distance_score_deriv = 0;
	Real distance_score = std::max<Real>( 0, std::cos( ( distance - 3.7 ) / 1.4 ) );
	if ( distance > 3.7 && distance < 6 ) {
		if ( distance_score > 0 ) {
			distance_score_deriv = - std::sin( ( distance - 3.7 ) / 1.4 ) / 1.4;
		}
	}

	Real dscore_ddistance = -1 * cylinder_score * plane_score *  distance_score_deriv;


	dscore_ddistance /= ring_atoms.size();

	Real scorefxn_weight = weights[arg_cation_pi];

	if ( dscore_ddistance != 0 ) {

		// All the derivative code is written assuming that ARG is res1
		// If it's not, we need to multiply by -1
		Real order_factor = arg_is_1 ? 1 : -1;

		Size arg_cz_ind = arg.atom_index( "CZ" );

		for ( std::string const & atom : ring_atoms ) {
			Size atom_ind = pi.atom_index(atom);
			auto xyz = pi.xyz( atom_ind );

			auto f1 = arg_cz.cross( xyz );
			auto f2 = ( arg_cz - xyz );

			auto f1s = dscore_ddistance * scorefxn_weight * f1;
			auto f2s = dscore_ddistance * scorefxn_weight * f2;

			Size atom1 = arg_is_1 ? arg_cz_ind : atom_ind;
			Size atom2 = arg_is_1 ? atom_ind : arg_cz_ind;

			r1_atom_derivs[ atom1 ].f1() += order_factor * f1s;
			r1_atom_derivs[ atom1 ].f2() += order_factor * f2s;
			r2_atom_derivs[ atom2 ].f1() -= order_factor * f1s;
			r2_atom_derivs[ atom2 ].f2() -= order_factor * f2s;
		}
	}
}


bool
ArgCationPiEnergy::defines_intrares_energy( EnergyMap const & ) const {
	return false;
}

void
ArgCationPiEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {
}


Distance
ArgCationPiEnergy::atomic_interaction_cutoff() const { return 6.0; }

core::Size
ArgCationPiEnergy::version() const { return 1; }


} // methods
} // scoring
} // core
